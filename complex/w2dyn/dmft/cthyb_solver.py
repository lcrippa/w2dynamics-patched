"""Module providing the high-level interface to the w2dynamics Fortran
CT-HYB impurity solver.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import time
from warnings import warn
import numpy as np
import scipy.optimize as opt
import h5py as h5
import sys
import ctypes as _c
# import copy

import w2dyn.auxiliaries.bindings as forlib
import w2dyn.auxiliaries.transform as tf
import w2dyn.auxiliaries.postprocessing as postproc
import w2dyn.auxiliaries.compound_index as ci
import w2dyn.auxiliaries.statistics as statistics
from w2dyn.auxiliaries.statistics import (DistributedSample,
                                          DistributedJackknife)
import w2dyn.dmft.dynamicalU as dynamicalU
from w2dyn.dmft.worm import (get_sector_index,
                             MAX_SECTOR,
                             SECTOR_Z,
                             SECTOR_CUSTOM)
from w2dyn.dmft.meta import QUANTITIES as qmetainfo
from w2dyn.dmft.impurity import (ImpuritySolver,
                                 StatisticalImpurityResult)
import w2dyn.dmft.orbspin as orbspin


class CtHybConfig:
    def __init__(self, taus, orbs, spins, cas, hashybs,
                 outer_sst):
        self.taus = taus
        self.orbs = orbs
        self.spins = spins
        self.cas = cas
        self.hashybs = hashybs
        self.outer_sst = outer_sst

    def __str__(self):
        return ("Monte Carlo configuration:\n" +
                "outer sst: {}\n".format(self.outer_sst) +
                "operators: {}".format(np.vstack((self.taus,
                                                  self.orbs,
                                                  self.spins,
                                                  self.cas,
                                                  self.hashybs))))

    @classmethod
    def get_from_ctqmc(cls):
        return cls(*forlib.get_qmc_config())

    def set_to_ctqmc(self):
        forlib.set_qmc_config(self.outer_sst,
                              self.taus, self.orbs, self.spins, self.cas,
                              self.hashybs)

    @classmethod
    def load_from_file(cls, file):
        with np.load(file) as f:
            return cls(f["taus"],
                       f["orbs"],
                       f["spins"],
                       f["cas"],
                       f["hashybs"],
                       f["outer_sst"][0])

    def save_to_file(self, file):
        np.savez(file,
                 taus=self.taus,
                 orbs=self.orbs,
                 spins=self.spins,
                 cas=self.cas,
                 hashybs=self.hashybs,
                 outer_sst=np.array([self.outer_sst]))


class CtHybSolver(ImpuritySolver):
    """Wrapper for Fortran CT-HYB(M-V) solver.

    This class provides access to the CT-QMC solver in the
    hybridisation expansion implemented using the matrix-vector
    technique by calling the native functions in the compiled Fortran
    library using the forlib module.

    Note that by being a Fortran module, it effectively makes this
    class a singleton, so you cannot run multiple instances of
    CtHybSolver on a single core in parallel.
    """
    def __init__(self, config, seed, Uw, Uw_Mat, epsn, interactive=True, mpi_comm=None):
        super(CtHybSolver, self).__init__(config)
        self.mpi_comm = mpi_comm
        self.mpi_rank = 0
        if mpi_comm is not None:
            self.mpi_rank = mpi_comm.Get_rank()

        self.seed = seed + self.mpi_rank
        self.param_string = self.config_to_pstring(config["QMC"])
        self.g_diagonal_only = (config["QMC"]["offdiag"] == 0)
        self.Uw = Uw
        self.Uw_Mat = Uw_Mat
        self.epsn = epsn

        # self.worm_config stores the worm configuration between the
        # processing of worm-related configuration during setup for a
        # single ctqmc_measure call and the extraction of worm-related
        # results after that call
        self.worm_config = None

        forlib.init_qmc_rng(seed)
        forlib.set_fancy_qmc_progress_bar(interactive)

        # START DYNAMICAL U
        if self.Uw == 1:
            Nd = config["Atoms"]["1"]["Nd"]
            beta = config["General"]["beta"]
            Nftau = config["QMC"]["Nftau"]
            umatrixdummy = np.zeros((Nd*2,Nd*2,Nd*2,Nd*2))
            self._Retarded2Shifts = dynamicalU.Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials(beta, umatrixdummy, self.Uw_Mat)
            screeningdummy = np.zeros((Nd,2,Nd,2,Nftau))
            self.screening = self._Retarded2Shifts.create_tau_dependent_UV_screening_function(Nftau, beta, screeningdummy)
            print("         ****************************")
            print("         ******* U(w) Message *******")
            print("         ****************************")
            print(" ")

    def prepare_quantum_numbers(self):
        """Calculates good quantum numbers based on self.problem.muimp
        and self.problem.interaction.auto_quantum_numbers and prints
        and returns the automatic quantum numbers if the default entry
        'auto' is contained in the interaction's quantum_numbers,
        otherwise prints a warning containing both automatic and
        manual quantum numbers and returns the manually specified
        quantum numbers."""

        # begin with qns permitted by interaction
        autoqns = self.problem.interaction.auto_quantum_numbers

        # threshold is relative to magnitude of largest entry in local
        # Hamiltonian
        threshold = max(
            np.amax(np.abs(self.problem.muimp)),
            np.amax(np.abs(self.problem.interaction.u_matrix))
        ) * abs(self.config["QMC"]["EPSBLOCK"])
        maxignored = 0.0
        elmsBelow = False

        if not self.g_diagonal_only:
            # check for violations by quadratic part of Hamiltonian
            spinodblocks = [np.abs(self.problem.muimp[:, s1, :, s2])
                            for s1 in range(self.problem.muimp.shape[1])
                            for s2 in range(self.problem.muimp.shape[3])
                            if s1 != s2]
            if any(np.any(block != 0) for block in spinodblocks):
                if any(np.any(block > threshold) for block in spinodblocks):
                    # spin-offdiagonals clear Azt, Szt, and (not necessarily)
                    # Jzt from good quantum numbers
                    if (any(qn in autoqns for qn in ("Azt", "Szt", "Jzt"))
                       and "All" not in autoqns):
                        # auto-partitioning might be worthwhile now
                        autoqns = autoqns + ("All", )
                    autoqns = tuple(qn for qn in autoqns
                                    if qn not in ("Azt", "Szt", "Jzt"))
                else:
                    maxignored = max(maxignored,
                                     *(np.amax(block)
                                       for block in spinodblocks))
                    elmsBelow = True
            orbodblocks = [np.abs(self.problem.muimp[o1, :, o2, :])
                           for o1 in range(self.problem.muimp.shape[0])
                           for o2 in range(self.problem.muimp.shape[2])
                           if o1 != o2]
            if any(np.any(block != 0) for block in orbodblocks):
                if any(np.any(block > threshold) for block in orbodblocks):
                    # orbital-offdiagonals clear Azt, Qzt, Lzt and (not
                    # necessarily) Jzt from good quantum numbers
                    if (any(qn in autoqns for qn in ("Azt", "Qzt", "Lzt", "Jzt"))
                       and "All" not in autoqns):
                        # auto-partitioning might be worthwhile now
                        autoqns = autoqns + ("All", )
                    autoqns = tuple(qn for qn in autoqns
                                    if qn not in ("Azt", "Qzt", "Lzt", "Jzt"))
                else:
                    maxignored = max(maxignored,
                                     *(np.amax(block)
                                       for block in orbodblocks))
                    elmsBelow = True

        if "auto" in self.problem.interaction.quantum_numbers:
            qns = " ".join(autoqns)
            if self.mpi_rank == 0:
                print("Setting automatic quantum numbers '{}'".format(qns),
                      file=sys.stderr)
                if elmsBelow:
                    warn("Off-diagonal elements of the local Hamiltonian below "
                         "EPSBLOCK (max. entry {:e}) were ignored for "
                         "automatic choice of quantum numbers"
                         .format(maxignored))
        else:
            if self.mpi_rank == 0:
                print("Overriding automatic quantum numbers '{}' with '{}'"
                      .format(
                          " ".join(autoqns),
                          " ".join(self.problem.interaction.quantum_numbers)
                      ), file=sys.stderr)
            if any(qn not in autoqns
                   for qn in self.problem.interaction.quantum_numbers
                   if qn != "All"):
                warn("Manually set quantum number(s) '{}' might be overly "
                     "restrictive;\n note that non-conforming elements of the "
                     "Hamiltonian will be IGNORED if present!"
                     .format(" ".join(
                         qn
                         for qn in self.problem.interaction.quantum_numbers
                         if qn not in autoqns and qn != "All"
                     )))
            qns = " ".join(self.problem.interaction.quantum_numbers)
        return qns

    def set_problem(self, problem, compute_fourpoint=0):
        # problem
        self.problem = problem
        problem_params = self.config_to_pstring({
            "Hamiltonian": problem.interaction.name,
            "beta": problem.beta,
            "Nd": problem.norbitals,
            "ParaMag": int(problem.paramag),
            "QuantumNumbers": self.prepare_quantum_numbers(),
            "Phonon": int(problem.use_phonons),
            "g_phonon": " ".join(problem.phonon_g),
            "omega0": " ".join(problem.phonon_omega0),
            "Screening": int(problem.use_screening),
            "Uw": self.Uw,
            })
        symm_params = self.symm_to_pstring(problem.symmetry_moves)
        all_params = self.param_string + problem_params + symm_params
        forlib.init_qmc_parameters(all_params)

        # remove hybridization and one-particle Hamiltonian off-diagonals in diagonal calculation
        if self.g_diagonal_only:
            self.ftau = orbspin.promote_diagonal(orbspin.extract_diagonal(self.problem.ftau))
            self.muimp = orbspin.promote_diagonal(orbspin.extract_diagonal(self.problem.muimp))
        else:
            self.ftau = self.problem.ftau
            self.muimp = self.problem.muimp

        if self.config["QMC"]["complexHyb"] == 0:
            if not np.allclose(np.imag(self.ftau), 0,
                               atol=self.config["QMC"]["real_ftau_tolerance"]):
                raise ValueError(
                    "The imaginary time hybridization function has "
                    "imaginary parts greater than QMC.real_ftau_tolerance "
                    "but QMC.complexHyb is not set"
                )
            self.ftau = np.real(self.ftau)
        if self.config["QMC"]["complexLocTr"] == 0:
            if not np.allclose(np.imag(self.muimp), 0,
                               atol=self.config["QMC"]["real_ham_tolerance"]):
                raise ValueError(
                    "The local Hamiltonian has imaginary parts "
                    "greater than QMC.real_ham_tolerance "
                    "but QMC.complexLocTr is not set"
                )
            self.muimp = np.real(self.muimp)
            if not np.allclose(np.imag(self.problem.interaction.u_matrix), 0,
                               atol=self.config["QMC"]["real_ham_tolerance"]):
                raise ValueError(
                    "The local interaction has imaginary parts "
                    "greater than QMC.real_ham_tolerance "
                    "but QMC.complexLocTr is not set"
                )
            self.problem.interaction.u_matrix = np.real(
                self.problem.interaction.u_matrix
            )

        # move tau-axis to last
        self.ftau = self.ftau.transpose(1,2,3,4,0)

        # CT-HYB uses a different convention for muimp
        self.muimp = -self.muimp
        self.umatrix = self.problem.interaction.u_matrix.reshape(
                             self.problem.nflavours, self.problem.nflavours,
                             self.problem.nflavours, self.problem.nflavours)
        if self.Uw == 0:
            self.screening = self.problem.screening.transpose(1, 2, 3, 4, 0).real

        # START DYNAMICAL U
        if self.Uw == 1:
            #print " "
            #print "         ****************************"
            #print "         ****** U(w) is active ******"
            #print "         ****************************"
            ### ===> UPDATE This is now done in the init above.
            ### ===> UPDATE _Retarded2Shifts = dynamicalU.Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials(problem.beta, self.umatrix, self.Uw_Mat)

            ### ===> UPDATE The umatrix shifting is now performed in config.py.
            ### ===> UPDATE The umatrix shifting is now performed in config.py.
            ### ===> # This should be executed for all atoms in the first iteration
            ### ===> try:
            ### ===>    self.iteration_ID    # Not defined for first DMFT iteration
            ### ===> except:
            ### ===>    self.umatrix_original = copy.copy(self.umatrix)
            ### ===>    self.iteration_ID = 1
            ### ===> if (self.umatrix==self.umatrix_original).all():
            ### ===>    self.umatrix = _Retarded2Shifts.shift_instantaneous_densitydensity_potentials(self.umatrix)
            ### ===> UPDATE The umatrix shifting is now performed in config.py.
            ### ===> UPDATE The umatrix shifting is now performed in config.py.

            ### ===> UPDATE This is now done in the init above.
            ### ===> UPDATE self.screening = _Retarded2Shifts.create_tau_dependent_UV_screening_function(problem.nftau, problem.beta, self.screening)

            # Perform chemical potential shift in each DMFT iteration
            if problem.norbitals == 1 and self.epsn == 0:
                self.muimp = self._Retarded2Shifts.shift_chemical_potential_for_single_band_case(self.muimp, self.umatrix)
            #else:
                #print "  ==> U(w): No shift of mu for mu-search"
                #if problem.norbitals > 1 and self.epsn==0:
                #print "===>>> U(w): mu-Search (epsn > 0) is recommended for the chosen problem."
        # END DYNAMICAL U


    def solve(self, iter_no, mc_config_inout=None):
        self.print_solver_config()
        sys.stdout.flush()
        sys.stderr.flush()
        if self.mpi_comm is not None:
            self.mpi_comm.Barrier()

        # Call CT-QMC solver
        forlib.set_qmc_simid(self.mpi_rank)
        time_qmc = time.perf_counter()
        #print('self.ftau', self.ftau)
        #exit()
        forlib.init_qmc_solver(self.umatrix, self.ftau,
                               self.muimp, self.screening)

        self.prepare_accumulators()

        firstrun = True
        if (self.config["QMC"]["ReuseMCConfig"] != 0
                and type(mc_config_inout) is list
                and len(mc_config_inout) > 0):
            firstrun = False
            mc_config_inout.pop(0).set_to_ctqmc()

        if self.config["QMC"]["segment"] == 0:
            taudiffmax = (self.config["QMC"]["TaudiffMax"]
                          if self.config["QMC"]["TaudiffMax"] > 0.0
                          else self.calibrate_taudiff_max())
            self.set_taudiffmax(taudiffmax)

        forlib.run_ctqmc_warmup(self.config["QMC"]["Nwarmups"]
                                if firstrun or self.config["QMC"]["Nwarmups2Plus"] < 0
                                else self.config["QMC"]["Nwarmups2Plus"])
        forlib.run_ctqmc_measurements(1, 1, self.config["QMC"]["MaxMeasurementTime"])
        time_qmc = time.perf_counter() - time_qmc
        if type(mc_config_inout) is list:
            mc_config_inout.append(CtHybConfig.get_from_ctqmc())

        aborted = forlib.get_qmc_aborted()
        if aborted == 2:  # no data, kill run
            raise KeyboardInterrupt("CTQMC was interrupted by signal")
        elif aborted == 1:
            self.abort = True
            warn("CT-QMC was aborted before all measurements were conducted.\n"
                 "Expect larger errorbars.", UserWarning, 2)

        # Get result set from the module variables.
        result = self.get_resultset(self.config["QMC"])
        result["time-qmc"] = time_qmc

        result.update(self.get_accumulator_results())

        rank_data_usable = (forlib.get_sector_meas_count(SECTOR_Z) > 0
                            and result["sign"] != 0.0)

        # Extract Green's function to use for Dyson equation
        giw = self.get_giw(self.config, self.problem, result)
        giw = DistributedSample([giw] if rank_data_usable else [],
                                self.mpi_comm)

        # Extract improved estimators
        gsigmaiw = self.get_gsigmaiw(self.config, self.problem, result)
        if gsigmaiw is not None:
            gsigmaiw = DistributedSample([gsigmaiw] if rank_data_usable else [],
                                         self.mpi_comm)

        # Extract moments of the self-energy from 1- and 2-particle reduced
        # density matrix
        try:
            rho1 = result["rho1"]
            rho2 = result["rho2"]
        except KeyError:
            smom = None
        else:
            umatrix_ctr = self.problem.interaction.u_matrix.reshape(
                                                [self.problem.nflavours] * 4)
            smom = postproc.get_siw_mom(umatrix_ctr, rho1, rho2)
            smom = DistributedSample([smom] if rank_data_usable else [],
                                     self.mpi_comm)

        result = {key: DistributedSample([value] if rank_data_usable else [],
                                         self.mpi_comm)
                  for (key, value) in result.items()}

        # Construct result object from result set and collect it from all nodes
        result = StatisticalImpurityResult(self.problem, giw, gsigmaiw, smom,
                                           self.g_diagonal_only, **result)

        # Cleanup
        forlib.dest_ctqmc()
        forlib.dest_qmc_paras()

        return result


    def solve_component(self,iter_no,isector,icomponent,mc_config_inout=[]):
        """ This solver returns the worm estimator for a given component.
           This results in several advantages:

              1. the balancing between worm space and partition function space
                 (i.e. eta-search) is done for specified component, and not in
                 an integrated way

              2. memory requirements are lowered for 2P quantities
                 (averaging can be done in python)

              3. (following 2.) no need for task-local files, postprocessing
                 formula to calculate H2->chi_conn can be implemented directly

           Convention of sampling spaces (Sector):

           Sector=1: Z (partition function, always sampled)
           Sector=2: G1 (one particle Green's function)
           Sector=3: GSigma (one particle improved estimator)
           Sector=4: G2 (two particle Green's function)
           Sector=5: H2 (two particle improved estimator)
           Sector=6: P2PH (single-freqeuncy 2P-GF in PH channel)
           Sector=7: P2PP (single-freqeuncy 2P-GF in PP channel)
           Sector=8: P3PH (two-freqeuncy 2P-GF in PH channel)
           Sector=9: P3PP (two-freqeuncy 2P-GF in PP channel)
           Sector=10: QQ (one-particle symmetric IE)
           Sector=11: QQQQ (two-particle symmetric IE)
           Sector=12: NQQdag
           Sector=13: QQdd
           Sector=14: Ucaca (P2 with U matrices)
           Sector=15: Uccaa (P2pp with U matrices)
           Sector=16: QUDdag
        """
        self.print_solver_config()
        sys.stdout.flush()
        sys.stderr.flush()
        if self.mpi_comm is not None:
            self.mpi_comm.Barrier()

        # Call CT-QMC solver
        forlib.set_qmc_simid(self.mpi_rank)
        time_qmc = time.perf_counter()
        forlib.init_qmc_solver(self.umatrix, self.ftau,
                               self.muimp, self.screening)
        forlib.set_wormeta(0.)

        self.prepare_accumulators()

        if type(mc_config_inout) is list and len(mc_config_inout) > 0:
            mc_config_inout[0].set_to_ctqmc()

        taudiffmax = (self.config["QMC"]["TaudiffMax"]
                      if self.config["QMC"]["TaudiffMax"] > 0.0
                      else self.calibrate_taudiff_max())
        self.set_taudiffmax(taudiffmax)

        # first warm-up
        if not mc_config_inout:
            forlib.run_ctqmc_warmup(self.config["QMC"]["Nwarmups"])
            mc_config_inout.append(CtHybConfig.get_from_ctqmc())

        if self.config["QMC"]["Nwarmups2Plus"] < 0:
            warn("Nwarmups2Plus not set;"
                 " assuming 0.1*Nwarmups for worm warmups and eta search steps")
            WormWarmups = self.config["QMC"]["Nwarmups"] // 10
        else:
            WormWarmups=self.config["QMC"]["Nwarmups2Plus"]

        only_zero_wormsteps = True
        if int(self.config['QMC']['WormSearchEta']) == 1:
            # Patrik's eta search: space step ratio indicates needed change of space eta ratio
            forlib.set_wormeta(np.double(self.config['QMC']['WormEta']),
                               get_sector_index(self.config["QMC"]))
            max_iter = 50
            autoraise_factor = 1.5
            n_local_extrema = 0
            nraise_local_extrema = 5
            eta_history = []

            for i_try in range(max_iter):
                steps_worm, steps_z = forlib.count_qmc_sector_steps(
                    isector,
                    icomponent,
                    WormWarmups,
                    WormWarmups
                )
                if not (steps_worm + steps_z == WormWarmups):
                    raise ValueError('{} steps in worm, {} steps in Z, '
                                     'some steps missing!'.format(steps_worm, steps_z))
                if steps_worm != 0:
                    only_zero_wormsteps = False
                res = abs(float(steps_worm - steps_z)) / float(steps_worm + steps_z)
                print(
                    'At eta={}: {} steps in worm, {} steps in Z, rel. diff. {}'.format(
                        forlib.get_wormeta(get_sector_index(self.config["QMC"])),
                        steps_worm, steps_z, res
                    )
                )

                # automatically raise step numbers if we are still
                # searching after having counted many local extrema in
                # eta
                eta_history.append(
                    forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                )
                if len(eta_history) >= 3:
                    if ((eta_history[-2] > eta_history[-3]
                         and eta_history[-2] > eta_history[-1])
                        or
                        (eta_history[-2] < eta_history[-3]
                         and eta_history[-2] < eta_history[-1])):
                        n_local_extrema += 1
                    if n_local_extrema >= nraise_local_extrema:
                        n_local_extrema = 0
                        WormWarmups = int(autoraise_factor * WormWarmups)

                if res > 0.1:
                    if i_try == max_iter - 1:
                        if only_zero_wormsteps:
                            print('No suitable value of eta found, setting it to 0.')
                            forlib.set_wormeta(0.)
                        else:
                            print('WARNING: Failed to reach eta search target accuracy.')
                        break
                    else:
                        forlib.set_wormeta(
                            forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                            * float(max(1, steps_z))
                            / float(max(1, steps_worm)),
                            get_sector_index(self.config["QMC"])
                        )
                else:
                    break
            eta_stat = DistributedSample(
                [forlib.get_wormeta(get_sector_index(self.config["QMC"]))],
                self.mpi_comm
            )
            eta_mean, eta_err = eta_stat.mean(), eta_stat.stderr()
            print('final eta on core: {}'.format(
                forlib.get_wormeta(get_sector_index(self.config["QMC"])))
            )
            print('final eta={}, stdev={}'.format(eta_mean, eta_err))

        elif int(self.config['QMC']['WormSearchEta']) == 2:
            # Josef's eta search: brentq root finding
            def eta_root(eta):
                '''
                Eta is determined by the root of this function.
                The ctqmc module function count_qmc_sector_steps is called
                to get the number of steps in Z and worm space
                for a given eta.
                '''
                forlib.set_wormeta(eta, get_sector_index(self.config["QMC"]))
                steps_worm, steps_z = forlib.count_qmc_sector_steps(
                    isector,
                    icomponent,
                    WormWarmups,
                    WormWarmups
                )
                if not (steps_worm + steps_z == WormWarmups):
                    raise ValueError('{} steps in worm, {} steps in Z, '
                                     'some steps missing!'.format(steps_worm, steps_z))
                res = float(steps_worm - steps_z) / float(steps_worm + steps_z)
                print('At eta={}: {} steps in worm, {} steps in Z -> {}'.format(eta, steps_worm, steps_z, res))
                return res

            bracket_min = 1e-4
            bracket_max = 1e2
            try:
                eta = opt.brentq(eta_root, bracket_min, bracket_max, xtol=0.001, rtol=0.001, maxiter=50, disp=True)
            except ValueError:
                try:  # for practical purposes, enlarging once should be sufficient.
                    print('Trying larger bracket for eta-search.')
                    eta = opt.brentq(eta_root, 0.01 * bracket_min, 100. * bracket_max, xtol=0.001, rtol=0.001, maxiter=50, disp=True)
                except ValueError:
                    print('No suitable value of eta found, setting it to 0.')
                    eta = 0.

            eta_stat = DistributedSample([eta], self.mpi_comm)
            eta_mean, eta_err = eta_stat.mean(), eta_stat.stderr()
            print('final eta on core: {}'.format(eta))
            print('final eta={}, stdev={}'.format(eta_mean, eta_err))

        else:
            eta_mean, eta_err = float(self.config['QMC']['WormEta']), 0.

        forlib.set_wormeta(eta_mean, get_sector_index(self.config["QMC"]))

        # skip unnecessary parts if worms cannot be inserted
        if eta_mean > 0.0:
            print('in python:')
            print('self.config["QMC"]["MaxMeasurementTime"] ', self.config["QMC"]["MaxMeasurementTime"] )
            forlib.run_ctqmc_worm_warmup(int(WormWarmups), isector, icomponent)
            forlib.run_ctqmc_measurements(isector, icomponent,
                                          self.config["QMC"]["MaxMeasurementTime"])

        time_qmc = time.perf_counter() - time_qmc

        aborted = forlib.get_qmc_aborted()
        if aborted == 2:  # no data, kill run
            raise KeyboardInterrupt("CTQMC was interrupted by signal")
        elif aborted == 1:
            self.abort = True
            warn("CT-QMC was aborted before all measurements were conducted.\n"
            "Expect larger errorbars.", UserWarning, 2)


        # partition space quantities without component tag
        result_gen = self.get_resultset(self.config["QMC"])
        result_gen.update(self.get_accumulator_results().items())

        rank_Z_usable = (forlib.get_sector_meas_count(SECTOR_Z) > 0
                         and result_gen["sign"] != 0.0)
        rank_worm_usable = (rank_Z_usable
                            and forlib.get_sector_meas_count(isector) > 0)

        result_gen = {key: DistributedSample([value] if rank_Z_usable else [],
                                             self.mpi_comm)
                      for key, value in result_gen.items()}

        result_gen = StatisticalImpurityResult(self.problem, None, None, None,
                                               self.g_diagonal_only,
                                               **result_gen)

        # partition space quantities with component tag + worm space
        # quantities (also with component tag)
        # FIXME: decide whether we really want to have the Z quants in
        # both (this allows checking data per component when they
        # would usually be averaged as it is already done in
        # solve_worm or overwritten as it would still be done in cthyb
        # and DMFT.py)
        result_comp = {key + "/{:05}".format(icomponent): value
                       for key, value in result_gen.other.items()}

        result_comp.update(
            {key + "/{:05}".format(icomponent):
             DistributedSample([value] if rank_worm_usable else [],
                               self.mpi_comm)
             for key, value
             in self.get_resultset_worm(self.config['QMC']).items()}
        )

        result_comp.update(
            {key + "/{:05}".format(icomponent):
             DistributedSample([value] if rank_worm_usable else [],
                               self.mpi_comm)
             for key, value in
             self.get_accumulator_results(worm=True).items()}
        )

        result_comp = StatisticalImpurityResult(self.problem, None, None, None,
                                                self.g_diagonal_only,
                                                **result_comp)

        # Cleanup
        forlib.dest_ctqmc()
        forlib.dest_qmc_paras()
        forlib.set_wormeta(0.)
        return result_gen, result_comp

    def solve_worm(self, iter_no, log_function=None):
        """
        Solve all non-vanishing components of G or GSigma or QQ by component sampling.
        Then convert the component list to standard format (matrix in orbitals).
        From this, compute the impurity Green's function and put it in a StatisticalImpurityResult.
        Furthermore, put all worm quantities to a separate StatisticalImpurityResult.
        Both of them are returned.
        :return: (StatisticalImpurityResult, StatisticalImpurityResult)
        """
        # get the components to be sampled
        conf_comp = self.config["QMC"]["WormComponents"]
        if conf_comp is not None: # use the ones specified in the parameters file
            components = map(int, conf_comp)
        elif self.config['QMC']['offdiag'] == 0: # generate all diagonal components
            components = [ci.GFComponent(bandspin=(iflav, iflav), n_ops=2, n_bands=self.problem.nflavours//2).index
                          for iflav in range(self.problem.nflavours)]
        else: # generate the full list
            components = 1 + np.arange(self.problem.nflavours**2) # one-based

        # introduce booleans for the measurement mode
        meas_g = self.config["QMC"]["WormMeasGiw"] != 0
        meas_gs = self.config["QMC"]["WormMeasGSigmaiw"] != 0
        meas_qq = self.config["QMC"]["WormMeasQQ"] != 0

        if meas_g and not meas_gs and not meas_qq:
            sector = 2
            # sector_name = 'giw-worm'
        elif not meas_g and meas_gs and not meas_qq:
            sector = 3
            # sector_name = 'gsigmaiw-worm'
        elif not meas_g and not meas_gs and meas_qq:
            sector = 10
            # sector_name = 'qqiw-worm'
        else:
            raise ValueError("Only one of (WormMeasGiw, "
                             "WormMeasGSigmaiw, WormMeasQQ) can be chosen")

        # perform a solve_component worm sampling run for each component.
        # in this run, the specified worm estimator is measured,
        # and also the densities in Z space, since we need them for moments etc.
        results_z = []
        result_list_worm = []
        conf_container = []
        worm_dict_complete = {}
        for component in components:
            log_function("sampling component {}".format(component))
            self.set_problem(self.problem)
            res_z, res_worm = self.solve_component(iter_no, sector, component,
                                                   conf_container)
            results_z.append((component, res_z))
            # FIXME: improve how we deal with missing worm component
            # results (when rank_worm_usable was False), but this will
            # require more changes below
            result_list_worm.append((component,
                                     {key: (value.local[0]
                                            if len(value.local) > 0
                                            else np.zeros_like(value.mean()))
                                      for key, value in res_worm.other.items()}))
            # merge worm results from all component runs (they are
            # distinguishable, since their names always contain the
            # component number, e.g. qqiw-worm/00001)
            worm_dict_complete.update(res_worm.other)

        # Z-sampling quantities, such as occ, are sampled together
        # with each component. To exploit all measurements, we collect
        # all samples together.
        def collect_components(qtty_name, result_list, mpi_comm):
            val_list = []
            for res in result_list:
                qtty = res[1].other[qtty_name]
                val_list.append(qtty)
            return statistics.join(*val_list)

        # average Z-sampled quantities from all worm component runs
        results_z = {key: collect_components(key, results_z, self.mpi_comm)
                     for key in results_z[0][1].other.keys()}

        # Extract moments of the self-energy from 1- and 2-particle reduced
        # density matrix
        umatrix_ctr = self.problem.interaction.u_matrix.reshape(
                                            [self.problem.nflavours] * 4)
        smom = [
            postproc.get_siw_mom(np.asarray(umatrix_ctr),
                                 np.asarray(rho1),
                                 np.asarray(rho2))
            for rho1, rho2 in zip(results_z['rho1'].local,
                                  results_z['rho2'].local)
        ]
        smom = DistributedSample(smom, self.mpi_comm)

        if meas_g:
            # convert component list of g to matrix
            giw = self.get_giw(self.config, self.problem, result_list_worm)
            gsigmaiw = None # not needed here
        if meas_gs:
            # convert component list of gsigma to matrix,
            # then apply IE equation to calculate giw
            gsigmaiw = self.get_gsigmaiw_worm(self.config, self.problem, result_list_worm)
            giw = self.giw_from_gsigmaiw_worm(self.config, self.problem, gsigmaiw)
        if meas_qq:
            # convert component list of qq to matrix,
            # then apply SIE equation to calculate giw
            gsigmaiw = None # not needed here
            qqiw = self.get_qqiw_worm(self.config, self.problem, result_list_worm)
            giw = self.giw_from_qqiw_worm(self.config, self.problem, qqiw, smom.local[0][0])

        giw = DistributedSample([giw], self.mpi_comm)
        if gsigmaiw is not None:
            gsigmaiw = DistributedSample([gsigmaiw], self.mpi_comm)

        return (StatisticalImpurityResult(self.problem, giw, gsigmaiw, smom, self.g_diagonal_only,
                                          **results_z),
                StatisticalImpurityResult(self.problem, giw, gsigmaiw, smom, self.g_diagonal_only,
                                          **worm_dict_complete))


        # -------- auxiliary helper methods ------------

    @classmethod
    def config_to_pstring(cls, config):
        return "".join("%s = %s #\n" % (key, value)
                       for key, value in sorted(config.items()))

    @classmethod
    def symm_to_pstring(cls, symmetries):
        s = "NSymMove = %d #\n" % len(symmetries)
        s +=  "".join("SymMove%02d = %s #\n" % (i+1, " ".join(map(str,vals+1)))
                      for i, vals in enumerate(symmetries))
        return s

    @classmethod
    def prepare_input_offdiag(cls, qtty):
        orbspin.warn_offdiagonal(qtty)
        #qtty = orbspin.extract_diagonal(qtty)
        if (np.abs(qtty.imag) > 1e-10).any():
            warn("Quantity has non-vanishing imaginary part", UserWarning, 2)
        qtty = qtty.real
        if qtty.ndim == 5:
            axes_list = [-4, -3, -2, -1] + list(range(qtty.ndim - 4))
            qtty = qtty.transpose(axes_list)
        return qtty

    @classmethod
    def prepare_input(cls, qtty):
        orbspin.warn_offdiagonal(qtty)
        qtty = orbspin.extract_diagonal(qtty)
        if (np.abs(qtty.imag) > 1e-10).any():
            warn("Quantity has non-vanishing imaginary part", UserWarning, 2)
        qtty = qtty.real
        if qtty.ndim > 2:
            axes_list = [-2, -1] + list(range(qtty.ndim - 2))
            qtty = qtty.transpose(axes_list)
        return qtty

    @classmethod
    def get_giw(cls, config, problem, result):
        #typ = "legendre"  # FIXME: QMC config does not contain FTtype
        typ = config["General"]["FTType"]
        if typ == "none":
            giw = -result["giw-full"]
        elif typ == "none_worm":
            giw = np.zeros_like(problem.g0inviw).transpose((1,2,3,4,0))
            for icomp, res in result:
                comp_ind = ci.GFComponent(index=icomp, n_ops=2,
                    n_bands=problem.nflavours//problem.nspins).bsbs()
                giw[comp_ind] = res.other['giw-worm/{:05}'.format(icomp)]
        elif typ == "plain":
            ntau = config["QMC"]["Ntau"]
            tf_matrix = tf.tau2mat(problem.beta, ntau, 'fermi',
                                   problem.niw)
            giw = tf.transform(tf_matrix, result["gtau-full"], None)
        elif typ == "legendre":
            nleg = config["QMC"]["NLegOrder"]
            tf_matrix = -tf.leg2mat(nleg, problem.niw)
            giw = tf.transform(tf_matrix, result["gleg-full"], None,
                               onmismatch=tf.Truncate.tofirst, warn=False)

            ### TODOcompl: the real version had this kind of symmetrization;
            # generalize it for complex F(tau) and t_ij
            #giw += giw.transpose(2, 3, 0, 1, 4)
            #giw /= 2.0

        else:
            raise RuntimeError("Invalid transform type: `%s'" % typ)

        # move freuency axis from back to front
        return np.ascontiguousarray(giw.transpose(4, 0, 1, 2, 3))

    @classmethod
    def get_gsigmaiw_worm(cls, config, problem, result):
        gsigma = np.zeros_like(problem.g0inviw).transpose((1,2,3,4,0))
        for icomp, res in result:
            comp_ind = ci.GFComponent(index=icomp, n_ops=2,
                n_bands=problem.nflavours//problem.nspins).bsbs()
            gsigma[comp_ind] = res.other['gsigmaiw-worm/{:05}'.format(icomp)]
        return gsigma.transpose((4,0,1,2,3))

    @classmethod
    def get_qqiw_worm(cls, config, problem, result):
        qq = np.zeros_like(problem.g0inviw).transpose((1, 2, 3, 4, 0))
        for icomp, res in result:
            comp_ind = ci.GFComponent(index=icomp, n_ops=2,
                    n_bands=problem.nflavours // problem.nspins).bsbs()
            qq[comp_ind] = res.other['qqiw-worm/{:05}'.format(icomp)]
        return qq.transpose((4, 0, 1, 2, 3))

    @classmethod
    def get_gsigmaiw(cls, config, problem, result):
        gsigmaiw = None
        # re-enable this when/if Z-space GSigma measurement is fixed
        # if config["QMC"]["MeasGSigmaiw"] == 1:
        #     gsigmaiw = result["gsigmaiw"]
        #     del result["gsigmaiw"] # have to delete it here, since gsigmaiw is passed explicitly to StatisticalImpurityResult
        return gsigmaiw

    @classmethod
    def giw_from_gsigmaiw_worm(cls, config, problem, gsigma):
        """
        Calculate impurity Green's function by improved estimator formula,
        Eq. (9) in PRB 100, 075119 (2019).
        """
        g0iw = orbspin.invert(problem.g0inviw)
        giw = g0iw + orbspin.multiply(g0iw, gsigma)
        return giw

    @classmethod
    def giw_from_qqiw_worm(cls, config, problem, qq, mom0):
        """
        Calculate impurity Green's function by symmetric improved estimator formula,
        Eq. (13) in PRB 100, 075119 (2019).
        """
        g0iw = orbspin.invert(problem.g0inviw)
        giw = g0iw + orbspin.multiply(g0iw, orbspin.multiply(qq + mom0, g0iw))
        return giw

    def calibrate_taudiff_max(self):
        phase1pairnum = 1000
        phase2taugrid = 100
        phase2pairnum = 1000
        nprogress = 100000
        forlib.reset_qmc_counters()
        acctaulist = np.zeros((phase1pairnum,), dtype=np.double)
        forlib.run_ctqmc_calibration(True,
                                     phase1pairnum, phase2taugrid, phase2pairnum,
                                     nprogress, acctaulist)
        acctaulist.sort()
        taudiffmax = DistributedSample([np.array((acctaulist[int(0.995 * phase1pairnum)],))],
                                       self.mpi_comm).mean()[0]
        self.set_taudiffmax(taudiffmax)
        forlib.reset_qmc_counters()
        accpair = np.zeros((self.problem.norbitals, 2, phase2taugrid), dtype=np.int64, order='F')
        forlib.run_ctqmc_calibration(False,
                                     phase1pairnum, phase2taugrid, phase2pairnum,
                                     nprogress, accpair)
        return self.estimate_taudiff_max(self.config,
                                         self.mpi_comm,
                                         accpair,
                                         taudiffmax,
                                         phase2taugrid)

    def set_taudiffmax(self, taudiffmax, log_function=None):
        if log_function is None:
            def log_function(*args, flush=True):
                if self.mpi_rank == 0:
                    print("CTQMC: ", *args, flush=flush)

        forlib.set_qmc_taudiffmax(taudiffmax)
        num_win = max(1, int(self.problem.beta / taudiffmax) - 1)
        log_function("Sliding window configuration: Using {num_win} window positions, "
                     "t_max = {tmax:.2f}, NWin = {nwin:.1f}, t_win = {twin:.2f}"
                     .format(num_win=num_win, tmax=taudiffmax, nwin=((num_win + 1)/2.0),
                             twin=self.problem.beta / ((num_win + 1)/2.0)))

    @staticmethod
    def estimate_taudiff_max(config, mpi_comm, accept_pair, taudiffmax, phase2taugrid):
        accept_pair = DistributedSample([accept_pair],
                                        mpi_comm).mean()

        def expdecay(p, taus, acc):
            return acc - p[0] * np.exp(p[1] * taus)

        estimates = np.empty((accept_pair.shape[0],
                              accept_pair.shape[1]),
                             dtype=np.double)
        for b in range(accept_pair.shape[0]):
            for s in range(accept_pair.shape[1]):
                fit = opt.leastsq(expdecay,
                                  (accept_pair[b, s, 0], -1),
                                  (tf.tau_bins(taudiffmax,
                                               phase2taugrid,
                                               "centre"),
                                   accept_pair[b, s]))
                estimates[b, s] = np.log(1 - 0.99) / fit[0][1]

        return np.amax(estimates)

    def print_solver_config(self, log_function=None):
        if log_function is None:
            def log_function(*args, flush=True):
                if self.mpi_rank == 0:
                    print("CTQMC: ", *args, flush=flush)

        if self.config["QMC"]["segment"] != 0:
            log_function("Using Segment solver")
        else:
            log_function("Using Matrix-Vector solver")
            if self.config["QMC"]["Eigenbasis"] == 0:
                raise ValueError("Krylov method requested in configuration is no longer supported")
            log_function("Using eigenbasis of impurity")
        if self.config["QMC"]["offdiag"] == 0:
            log_function("Using diagonal F(tau)")
        else:
            log_function("Using full F(tau) including offdiagonals")
            if self.config["QMC"]["full_offdiag"] != 0:
                log_function("Fast matrix updates disabled")
        if self.config["QMC"]["MeasSusz"] != 0:
            log_function("Measuring susceptibility")
        if self.config["QMC"]["MeasGSigmaiw"] != 0:
            log_function("WARNING: MeasGSigmaiw /= 0: Measuring GSigmaiw"
                         " in Z sampling is not currently supported!")

    def prepare_accumulators(self):
        # Get collection of measurement targets used by CTQMC to set our accumulators and samplers
        meastargets = forlib.get_global_meastargets_FIXME()

        self.accumulators = {}
        self.samplers = {}

        # Prepare accumulators and samplers for measurements
        self.accumulators['Z_G_tau'] = forlib.MeanAcc(
            self.problem.nflavours**2 * self.config["QMC"]["Ntau"],
            True
        )
        self.samplers['Z_G_tau'] = forlib.RHistSampler(
            self.problem.nflavours**2,
            self.problem.beta,
            -1,
            self.config["QMC"]["Ntau"],
            1
        )
        forlib.register_sampler_G(meastargets,
                                  self.samplers['Z_G_tau'],
                                  self.accumulators['Z_G_tau'])

        self.accumulators['Z_G_leg'] = forlib.MeanAcc(
            self.problem.nflavours**2 * self.config["QMC"]["NLegOrder"],
            True
        )
        self.samplers['Z_G_leg'] = forlib.LegendreSampler(
            self.problem.nflavours**2,
            self.problem.beta,
            -1,
            self.config["QMC"]["NLegOrder"],
            1
        )
        forlib.register_sampler_G(meastargets,
                                  self.samplers['Z_G_leg'],
                                  self.accumulators['Z_G_leg'])
        
        try:
            sparse_frequencies = list(np.loadtxt("sparse_frequencies.dat")) # loading IR sampling frequencies if given
        except:
            sparse_frequencies = []

        self.samplers['Z_G_mat'] = forlib.FourierSampler(
            self.problem.nflavours**2,
            self.problem.beta,
            -1,
            self.config["QMC"]["Niw"],
            1,
            use_nfft=True,
            sparse_frequencies=sparse_frequencies
        )
        self.accumulators['Z_G_mat'] = forlib.MeanAcc(
            self.samplers['Z_G_mat'].nacc * self.samplers['Z_G_mat'].ncomp,
            True
        )
        forlib.register_sampler_G(meastargets,
                                  self.samplers['Z_G_mat'],
                                  self.accumulators['Z_G_mat'])

        self.accumulators['Z_sign'] = forlib.MeanAcc(1, True)
        forlib.set_accumulator_sign(meastargets, self.accumulators['Z_sign'])


        if self.config["QMC"]["MeasG2iw"]:
            # sampler and accumulator for two-frequency GF 
            self.samplers['Z_G2iw_mat'] = forlib.Fourier2DSampler(
                self.problem.nflavours**2,
                self.problem.beta,
                [-1,-1],
                [self.config["QMC"]["Niw"], self.config["QMC"]["Niw"]],
                np.eye(2),
                use_nfft=True
            )
            self.accumulators['Z_G2iw_mat'] = forlib.MeanAcc(
                self.samplers['Z_G2iw_mat'].nacc * self.samplers['Z_G2iw_mat'].ncomp,
                True
            )
            forlib.register_sampler_G2iw(meastargets,
                                    self.samplers['Z_G2iw_mat'],
                                    self.accumulators['Z_G2iw_mat'])


        if self.config["QMC"]["FourPnt"]==8:
            try:
                file = h5.File("sampling_points.hdf5")
                sampling_points = file['sampling_points'][()].T
            except FileNotFoundError:
                raise FileNotFoundError("For the measurement of the general four-fermionic-frequency" \
                "two-particel Green's function the IR basis sampling points have to be provided in the 'sampling_points.hdf5' file.")

            del file

            self.accumulators['Z_G4iw_full'] = forlib.MeanAcc(
                self.problem.nflavours**4 * sampling_points.shape[-1],
                True
            )        
        
            sparse_frequencies = np.ascontiguousarray(np.reshape(sampling_points, 4*sampling_points.shape[-1]), np.int64)
            sparse_frequencies_arg = sparse_frequencies.ctypes.data_as(_c.POINTER(_c.c_int64))

            forlib.set_accumulator_G4iw_full(meastargets, self.accumulators['Z_G4iw_full'], sparse_frequencies_arg, np.shape(sparse_frequencies)[-1], \
                                             self.problem.norbitals, self.problem.beta)
            

        self.accumulators['Z_dm'] = forlib.MeanAcc(
            (2**self.problem.nflavours)**2,  # nstates x nstates
            True
        )
        forlib.set_accumulator_dm(meastargets, self.accumulators['Z_dm'])

        self.accumulators['expansion_order'] = forlib.MeanAcc(
            self.problem.nflavours * (self.config["QMC"]["MaxHisto"] + 1),
            True
        )
        forlib.set_accumulator_expord(meastargets, self.accumulators['expansion_order'])

        # XXX why is this not complex as the sign is?
        self.accumulators['Z_ntau_n0'] = forlib.MeanAcc(
            self.problem.nflavours**2 * (self.config["QMC"]["Ntau"] + 1),
            True
        )
        forlib.set_accumulator_ntaun0(meastargets, self.accumulators['Z_ntau_n0'])

        self.worm_config = None

        if (self.config["QMC"]["WormMeasGiw"]
            or self.config["QMC"]["WormMeasGtau"]):
            # FIXME: perhaps more sensible to group elements of stat,
            # conv, quant_mat, quant_tau into (ideally named) tuples
            # (one per result quantity), replace quant_x by type and
            # name entry in those tuples and maybe add normalization
            # factors
            self.worm_config = {
                'ntauops': (1, 1),
                'optypes': (forlib.OP_W.ANNH, forlib.OP_W.CREA),
                'stat': [],
                'conv': [],
                'quant_mat': [],
                'quant_tau': [],
            }

        if self.config["QMC"]["WormMeasGiw"]:
            self.worm_config['stat'].append(
                (forlib.STAT_FERMI,)
            )
            self.worm_config['conv'].append(
                # tau1 - tau2
                np.array(((1.0,),))
            )
            self.worm_config['quant_mat'].append(
                'giw-worm'
            )
            self.worm_config['quant_tau'].append(
                None
            )

        if self.config["QMC"]["WormMeasGtau"]:
            # FIXME: result has factor (-1) relative to Z G(tau),
            # change normalization? (requires ideally supporting
            # normalization factors in self.worm_config)
            self.worm_config['stat'].append(
                (forlib.STAT_FERMI,)
            )
            self.worm_config['conv'].append(
                # tau1 - tau2
                np.array(((1.0,),))
            )
            self.worm_config['quant_mat'].append(
                None
            )
            self.worm_config['quant_tau'].append(
                'gtau-worm'
            )

        if (self.config["QMC"]["WormMeasG4tau"]
            or self.config["QMC"]["WormMeasG4iw"]):
            if self.worm_config is not None:
                raise ValueError("Only one worm estimator "
                                 "can be enabled per run")
            self.worm_config = {
                'ntauops': (1, 1, 1, 1),
                'optypes': (forlib.OP_W.ANNH, forlib.OP_W.CREA,
                            forlib.OP_W.ANNH, forlib.OP_W.CREA),
                'stat': [],
                'conv': [],
                'quant_mat': [],
                'quant_tau': [],
            }

        if self.config["QMC"]["WormMeasG4iw"]:
            if self.config["QMC"]["G4ph"]:
                self.worm_config['stat'].append(
                    (forlib.STAT_FERMI,
                     forlib.STAT_FERMI,
                     forlib.STAT_BOSE)
                )
                self.worm_config['conv'].append(
                    # tau1 - tau2, tau3 - tau4, tau1 - tau4
                    np.array(((1.0, -1.0, 0.0),
                              (0.0, 0.0, 1.0),
                              (1.0, 0.0, 0.0)))
                    if self.config["QMC"]["WormPHConvention"] == 0
                    # tau1 - tau2, tau3 - tau4, tau2 - tau3
                    else np.array(((1.0, -1.0, 0.0),
                                   (0.0, 0.0, 1.0),
                                   (0.0, 1.0, -1.0))),
                )
                self.worm_config['quant_mat'].append(
                    'g4iw-worm'
                )
                self.worm_config['quant_tau'].append(
                    None
                )

            if self.config["QMC"]["G4pp"]:
                self.worm_config['stat'].append(
                    (forlib.STAT_FERMI,
                     forlib.STAT_FERMI,
                     forlib.STAT_BOSE)
                )
                self.worm_config['conv'].append(
                    # tau1 - tau3, tau2 - tau4, tau3 - tau2
                    np.array(((1.0, 0.0, -1.0),
                              (0.0, 1.0, 0.0),
                              (0.0, -1.0, 1.0)))
                )
                self.worm_config['quant_mat'].append(
                    'g4iwpp-worm'
                )
                self.worm_config['quant_tau'].append(
                    None
                )

        if self.config["QMC"]["WormMeasG4tau"]:
            self.worm_config['stat'].append(
                (forlib.STAT_FERMI,
                 forlib.STAT_FERMI,
                 forlib.STAT_BOSE)
            )
            self.worm_config['conv'].append(
                # tau1 - tau2, tau3 - tau4, tau1 - tau4
                np.array(((1.0, -1.0, 0.0),
                          (0.0, 0.0, 1.0),
                          (1.0, 0.0, 0.0)))
            )
            self.worm_config['quant_mat'].append(
                None
            )
            self.worm_config['quant_tau'].append(
                'g4tau-worm'
            )

        if self.config["QMC"]["WormMeasP2iwPH"]:
            self.samplers['Worm_P2ph_mat'] = forlib.FourierSampler(
                1,
                self.problem.beta,
                forlib.STAT_BOSE,
                self.config["QMC"]["N2iwb"],
                1
            )
            self.accumulators['Worm_P2ph_mat'] = forlib.MeanAcc(
                self.samplers['Worm_P2ph_mat'].nacc * self.samplers['Worm_P2ph_mat'].ncomp,
                True
            )
            forlib.register_sampler_wormP2ph(meastargets,
                                             self.samplers['Worm_P2ph_mat'],
                                             self.accumulators['Worm_P2ph_mat'])

        if self.config["QMC"]["WormMeasP2tauPH"]:
            self.samplers['Worm_P2ph_tau'] = forlib.RHistSampler(
                1,
                self.problem.beta,
                forlib.STAT_BOSE,
                self.config["QMC"]["Ntau"],
                1
            )
            self.accumulators['Worm_P2ph_tau'] = forlib.MeanAcc(
                self.samplers['Worm_P2ph_tau'].nacc * self.samplers['Worm_P2ph_tau'].ncomp,
                True
            )
            forlib.register_sampler_wormP2ph(meastargets,
                                             self.samplers['Worm_P2ph_tau'],
                                             self.accumulators['Worm_P2ph_tau'])

        if self.config["QMC"]["WormMeasP2iwPP"]:
            self.samplers['Worm_P2pp_mat'] = forlib.FourierSampler(
                1,
                self.problem.beta,
                forlib.STAT_BOSE,
                self.config["QMC"]["N2iwb"],
                1
            )
            self.accumulators['Worm_P2pp_mat'] = forlib.MeanAcc(
                self.samplers['Worm_P2pp_mat'].nacc * self.samplers['Worm_P2pp_mat'].ncomp,
                True
            )
            forlib.register_sampler_wormP2pp(meastargets,
                                             self.samplers['Worm_P2pp_mat'],
                                             self.accumulators['Worm_P2pp_mat'])

        if self.config["QMC"]["WormMeasP2tauPP"]:
            self.samplers['Worm_P2pp_tau'] = forlib.RHistSampler(
                1,
                self.problem.beta,
                forlib.STAT_BOSE,
                self.config["QMC"]["Ntau"],
                1
            )
            self.accumulators['Worm_P2pp_tau'] = forlib.MeanAcc(
                self.samplers['Worm_P2pp_tau'].nacc * self.samplers['Worm_P2pp_tau'].ncomp,
                True
            )
            forlib.register_sampler_wormP2pp(meastargets,
                                             self.samplers['Worm_P2pp_tau'],
                                             self.accumulators['Worm_P2pp_tau'])

        if (self.config["QMC"]["WormMeasP3iwPH"]
            or self.config["QMC"]["WormMeasP3tauPH"]):
            if self.worm_config is not None:
                raise ValueError("Only one worm estimator "
                                 "can be enabled per run")
            self.worm_config = {
                'ntauops': (1, 1, 2),
                'optypes': (forlib.OP_W.CREA, forlib.OP_W.ANNH,
                            forlib.OP_W.CREA, forlib.OP_W.ANNH),
                'stat': [],
                'conv': [],
                'quant_mat': [],
                'quant_tau': [],
            }

        if self.config["QMC"]["WormMeasP3iwPH"]:
            self.worm_config['stat'].append(
                (forlib.STAT_FERMI, forlib.STAT_BOSE)
            )
            self.worm_config['conv'].append(
                np.array(((1.0, -1.0),
                          (0.0, 1.0)))
            )
            self.worm_config['quant_mat'].append(
                'p3iw-worm'
            )
            self.worm_config['quant_tau'].append(
                None
            )

        if self.config["QMC"]["WormMeasP3tauPH"]:
            self.worm_config['stat'].append(
                (forlib.STAT_FERMI, forlib.STAT_BOSE)
            )
            self.worm_config['conv'].append(
                np.array(((1.0, -1.0),
                          (0.0, 1.0)))
            )
            self.worm_config['quant_mat'].append(
                None
            )
            self.worm_config['quant_tau'].append(
                'p3phtau-worm'
            )

        if (self.config["QMC"]["WormMeasP3iwPP"]
            or self.config["QMC"]["WormMeasP3tauPP"]):
            if self.worm_config is not None:
                raise ValueError("Only one worm estimator "
                                 "can be enabled per run")
            self.worm_config = {
                'ntauops': (1, 1, 2),
                'optypes': (forlib.OP_W.ANNH, forlib.OP_W.ANNH,
                            forlib.OP_W.CREA, forlib.OP_W.CREA),
                'stat': [],
                'conv': [],
                'quant_mat': [],
                'quant_tau': [],
            }

        if self.config["QMC"]["WormMeasP3iwPP"]:
            self.worm_config['stat'].append(
                (forlib.STAT_FERMI, forlib.STAT_BOSE)
            )
            self.worm_config['conv'].append(
                np.array(((1.0, -1.0),
                          (0.0, 1.0)))
            )
            self.worm_config['quant_mat'].append(
                'p3iwpp-worm'
            )
            self.worm_config['quant_tau'].append(
                None
            )

        if self.config["QMC"]["WormMeasP3tauPP"]:
            self.worm_config['stat'].append(
                (forlib.STAT_FERMI, forlib.STAT_BOSE)
            )
            self.worm_config['conv'].append(
                np.array(((1.0, -1.0),
                          (0.0, 1.0)))
            )
            self.worm_config['quant_mat'].append(
                None
            )
            self.worm_config['quant_tau'].append(
                'p3pptau-worm'
            )

        if self.config["QMC"]["WormMeasQQ"]:
            self.samplers['Worm_QQ_mat'] = forlib.FourierSampler(
                1,
                self.problem.beta,
                forlib.STAT_FERMI,
                self.config["QMC"]["Niw"],
                1
            )
            self.accumulators['Worm_QQ_mat'] = forlib.MeanAcc(
                self.samplers['Worm_QQ_mat'].nacc * self.samplers['Worm_QQ_mat'].ncomp,
                True
            )
            forlib.register_sampler_wormQQ(meastargets,
                                           self.samplers['Worm_QQ_mat'],
                                           self.accumulators['Worm_QQ_mat'])

        if self.config["QMC"]["WormMeasQQtau"]:
            self.samplers['Worm_QQ_tau'] = forlib.RHistSampler(
                1,
                self.problem.beta,
                forlib.STAT_FERMI,
                self.config["QMC"]["Ntau"],
                1
            )
            self.accumulators['Worm_QQ_tau'] = forlib.MeanAcc(
                self.samplers['Worm_QQ_tau'].nacc * self.samplers['Worm_QQ_tau'].ncomp,
                True
            )
            forlib.register_sampler_wormQQ(meastargets,
                                           self.samplers['Worm_QQ_tau'],
                                           self.accumulators['Worm_QQ_tau'])

        if self.config["QMC"]["WormMeasRaman"]:
            self.samplers['Worm_Raman_mat'] = forlib.Fourier2DSampler(
                1,
                self.problem.beta,
                [forlib.STAT_BOSE, forlib.STAT_BOSE],
                [self.config["QMC"]["N3iwb"], self.config["QMC"]["N3iwb"]],
                np.eye(2)
            )
            self.accumulators['Worm_Raman_mat'] = forlib.MeanAcc(
                self.samplers['Worm_Raman_mat'].nacc * self.samplers['Worm_Raman_mat'].ncomp,
                True
            )
            forlib.register_sampler2d_wormRaman(meastargets,
                                                self.samplers['Worm_Raman_mat'],
                                                self.accumulators['Worm_Raman_mat'])

        if (self.config["QMC"]["WormMeasCustomMat"]
            or self.config["QMC"]["WormMeasCustomTau"]
            or self.config["QMC"]["WormSampleCustom"]):

            if self.worm_config is not None:
                raise ValueError("Only one worm estimator "
                                 "can be enabled per run")

            ntauops = [int(i)
                       for i
                       in self.config["QMC"]["WormCustomNOpsTau"]]
            optypes = [int(i)
                       for i
                       in self.config["QMC"]["WormCustomOptypes"]]
            self.worm_config = {
                'ntauops': ntauops,
                'optypes': optypes,
            }
        elif self.worm_config is not None:
            ntauops = self.worm_config['ntauops']
            optypes = self.worm_config['optypes']

        if self.worm_config is not None:
            nwormtaus = len(ntauops)
            optypegroups = []
            if (sum(ntauops) != len(optypes)):
                raise ValueError("Custom worm: sum(NOpsTau) != len(Optypes)")
            for nops in ntauops:
                optypegroups.append(optypes[:nops])
                optypes = optypes[nops:]
            optypes = [i for g in optypegroups for i in g]
            forlib.set_custom_worm(
                forlib.get_wormstate(),
                optypegroups
            )
            opdesc = " ".join(
                ["".join(["c" if x == forlib.OP_W.ANNH else "c^\\dag"
                          for x in reversed(eqtauprod)])
                 + "(tau_{})".format(i + 1)
                 for i, eqtauprod in enumerate(optypegroups)]
            )


            if 'stat' in self.worm_config:
                stat = self.worm_config["stat"]
            else:
                stat = self.config["QMC"]["WormCustomStat"]
                if any(s not in ('fermionic', 'bosonic') for s in stat):
                    raise ValueError("Invalid entry in WormCustomStat")
                if len(stat) == 0 or len(stat) % (nwormtaus - 1) != 0:
                    raise ValueError("Invalid length of WormCustomStat")
                stat = [
                    [forlib.STAT_FERMI if x == 'fermionic' else forlib.STAT_BOSE
                     for x in stat[i * (nwormtaus - 1)
                                   : (i + 1) * (nwormtaus - 1)]]
                    for i in range(len(stat) // (nwormtaus - 1))
                ]
                self.worm_config['stat'] = stat

            if 'conv' in self.worm_config:
                conv = self.worm_config["conv"]
            elif self.config["QMC"]["WormCustomConvention"] is None:
                conv = [np.eye(nwormtaus - 1, dtype=np.double) for s in stat]
                self.worm_config['conv'] = conv
            else:
                conv = self.config["QMC"]["WormCustomConvention"]
                if (len(conv) == 0 or len(conv) % (nwormtaus - 1)**2 != 0):
                    raise ValueError("Invalid length of WormCustomConvention")
                conv = [
                    np.reshape(
                        [float(x)
                         for x in conv[i * (nwormtaus - 1)**2
                                       : (i + 1) * (nwormtaus - 1)**2]],
                        (nwormtaus - 1, nwormtaus - 1)
                    )
                    for i in range(len(conv) // (nwormtaus - 1)**2)
                ]
                self.worm_config['conv'] = conv

        if self.config["QMC"]["WormMeasCustomMat"]:
            if self.config["QMC"]["WormCustomMatQuant"] is None:
                if len(stat) == 1:
                    self.worm_config["quant_mat"] = ["custom-iw-worm"]
                else:
                    self.worm_config["quant_mat"] = [
                        "custom-iw-worm-conv{}".format(i+1)
                        for i in range(len(stat))
                    ]
            else:
                self.worm_config["quant_mat"] = [
                    str(x) for x in self.config["QMC"]["WormCustomMatQuant"]
                ]
        elif (self.worm_config is not None
              and 'quant_mat' not in self.worm_config):
            self.worm_config["quant_mat"] = []

        if self.config["QMC"]["WormMeasCustomTau"]:
            if self.config["QMC"]["WormCustomTauQuant"] is None:
                if len(stat) == 1:
                    self.worm_config["quant_tau"] = ["custom-tau-worm"]
                else:
                    self.worm_config["quant_tau"] = [
                        "custom-tau-worm-conv{}".format(i+1)
                        for i in range(len(stat))
                    ]
            else:
                self.worm_config["quant_tau"] = [
                    str(x) for x in self.config["QMC"]["WormCustomTauQuant"]
                ]
        elif (self.worm_config is not None
              and 'quant_tau' not in self.worm_config):
            self.worm_config["quant_tau"] = []

        if (self.worm_config is not None
            and any(q is not None for q in self.worm_config["quant_mat"])):
            for i, (stat, conv, qname) in enumerate(
                    zip(self.worm_config['stat'],
                        self.worm_config['conv'],
                        self.worm_config['quant_mat'])
            ):
                if qname is None:
                    continue
                key = 'Worm_Custom_mat{}'.format(i)
                # If a description already exists, use it as prefix
                if qname in qmetainfo:
                    desc = (qmetainfo[qname]['desc'] + ';\n'
                            '(automatically generated description:)\n')
                else:
                    desc = ''
                # FIXME: clearer description, esp. convention
                desc += (
                    "Worm estimator {} on Matsubara frequencies, "
                    "convention transformation {}"
                    .format(opdesc, conv)
                )

                if nwormtaus == 2:
                    qmetainfo[qname] = {
                        "axes": [
                            "ineq",
                            "iw" if stat[0] == forlib.STAT_FERMI else "iwb-p2"
                        ],
                        "desc": desc
                    }
                    self.samplers[key] = forlib.FourierSampler(
                        1,
                        self.problem.beta,
                        *stat,
                        self.config["QMC"]["Niw"]
                        if stat[0] == forlib.STAT_FERMI
                        # FIXME: consistent meaning of num. bosonic freq.
                        else (self.config["QMC"]["N2iwb"] + 1),
                        conv[0, 0]
                    )
                elif nwormtaus == 3:
                    qmetainfo[qname] = {
                        "axes": [
                            "ineq",
                            "iwf-p3" if stat[0] == forlib.STAT_FERMI else "iwb-p3",
                            "iwf-p3" if stat[1] == forlib.STAT_FERMI else "iwb-p3",
                        ],
                        "desc": desc
                    }
                    self.samplers[key] = forlib.Fourier2DSampler(
                        1,
                        self.problem.beta,
                        stat,
                        [self.config["QMC"]["N3iwf"]
                         if s == forlib.STAT_FERMI
                         # FIXME: consistent meaning of num. bosonic freq.
                         else (self.config["QMC"]["N3iwb"] + 1)
                         for s in stat],
                        conv
                    )
                elif nwormtaus == 4:
                    qmetainfo[qname] = {
                        "axes": [
                            "ineq",
                            "iwf-g4" if stat[0] == forlib.STAT_FERMI else "iwb-g4",
                            "iwf-g4" if stat[1] == forlib.STAT_FERMI else "iwb-g4",
                            "iwf-g4" if stat[2] == forlib.STAT_FERMI else "iwb-g4",
                        ],
                        "desc": desc
                    }
                    self.samplers[key] = forlib.Fourier3DSampler(
                        1,
                        self.problem.beta,
                        stat,
                        [self.config["QMC"]["N4iwf"]
                         if s == forlib.STAT_FERMI
                         # FIXME: consistent meaning of num. bosonic freq.
                         else (self.config["QMC"]["N4iwb"] + 1)
                         for s in stat],
                        conv
                    )
                else:
                    raise NotImplementedError("Measurement on Matsubara "
                                              "frequencies not supported for "
                                              "this number of custom worm taus")

                self.accumulators[key] = forlib.MeanAcc(
                    self.samplers[key].nacc
                    * self.samplers[key].ncomp,
                    True
                )
                forlib.register_sampler_wormCustom(
                    meastargets,
                    self.samplers[key],
                    self.accumulators[key],
                    nwormtaus
                )


        if (self.worm_config is not None
            and any(q is not None for q in self.worm_config["quant_tau"])):
            for i, (stat, conv, qname) in enumerate(
                    zip(self.worm_config['stat'],
                        self.worm_config['conv'],
                        self.worm_config['quant_tau'])
            ):
                if qname is None:
                    continue
                key = 'Worm_Custom_tau{}'.format(i)
                # If a description already exists, use it as prefix
                if qname in qmetainfo:
                    desc = (qmetainfo[qname]['desc'] + ';\n'
                            '(automatically generated description:)\n')
                else:
                    desc = ''
                # FIXME: clearer description, esp. convention
                desc += (
                    "Worm estimator {} on tau bins, "
                    "convention transformation {}"
                    .format(opdesc, conv)
                )

                if nwormtaus == 2:
                    qmetainfo[qname] = {
                        "axes": [
                            "ineq",
                            "taubin"
                        ],
                        "desc": desc
                    }
                    self.samplers[key] = forlib.RHistSampler(
                        1,
                        self.problem.beta,
                        *stat,
                        self.config["QMC"]["Ntau"],
                        conv[0, 0]
                    )
                elif nwormtaus == 3:
                    qmetainfo[qname] = {
                        "axes": [
                            "ineq",
                            "taubin",
                            "taubin",
                        ],
                        "desc": desc
                    }
                    self.samplers[key] = forlib.RHist2DSampler(
                        1,
                        self.problem.beta,
                        stat,
                        [self.config["QMC"]["Ntau"]
                         for s in stat],
                        conv
                    )
                elif nwormtaus == 4:
                    qmetainfo[qname] = {
                        "axes": [
                            "ineq",
                            "tau-g4",
                            "tau-g4",
                            "tau-g4",
                        ],
                        "desc": desc
                    }
                    self.samplers[key] = forlib.RHist3DSampler(
                        1,
                        self.problem.beta,
                        stat,
                        [self.config["QMC"]["N4tau"]
                         for s in stat],
                        conv
                    )
                else:
                    raise NotImplementedError("Measurement into tau bins "
                                              "not supported for "
                                              "this number of custom worm taus")
                self.accumulators[key] = forlib.MeanAcc(
                    self.samplers[key].nacc
                    * self.samplers[key].ncomp,
                    True
                )
                forlib.register_sampler_wormCustom(
                    meastargets,
                    self.samplers[key],
                    self.accumulators[key],
                    nwormtaus
                )

    def get_accumulator_results(self, worm=False):
        result = {}

        mean_sign = self.accumulators['Z_sign'].mean[0]

        if worm:

            if self.config["QMC"]["WormMeasP2iwPH"]:
                # FIXME: normalization (check/fix)
                result["p2iw-worm"] = (
                    forlib.get_sector_sample_count(get_sector_index(self.config["QMC"]))
                    / forlib.get_sector_sample_count(SECTOR_Z)
                    / self.problem.norbitals**4
                    / forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                    * self.samplers['Worm_P2ph_mat']
                    .postprocess(self.accumulators['Worm_P2ph_mat'].mean)[0, ...]
                ) / mean_sign

            if self.config["QMC"]["WormMeasP2tauPH"]:
                # FIXME: normalization (check/fix)
                result["p2tau-worm"] = (
                    forlib.get_sector_sample_count(get_sector_index(self.config["QMC"]))
                    / forlib.get_sector_sample_count(SECTOR_Z)
                    / self.problem.norbitals**4
                    / forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                    * self.samplers['Worm_P2ph_tau']
                    .postprocess(self.accumulators['Worm_P2ph_tau'].mean)[0, ...]
                ) / mean_sign

            if self.config["QMC"]["WormMeasP2iwPP"]:
                # FIXME: normalization (check/fix)
                result["p2iwpp-worm"] = (
                    forlib.get_sector_sample_count(get_sector_index(self.config["QMC"]))
                    / forlib.get_sector_sample_count(SECTOR_Z)
                    / self.problem.norbitals**4
                    / forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                    * self.samplers['Worm_P2pp_mat']
                    .postprocess(self.accumulators['Worm_P2pp_mat'].mean)[0, ...]
                ) / mean_sign

            if self.config["QMC"]["WormMeasP2tauPP"]:
                # FIXME: normalization (check/fix)
                result["p2taupp-worm"] = (
                    forlib.get_sector_sample_count(get_sector_index(self.config["QMC"]))
                    / forlib.get_sector_sample_count(SECTOR_Z)
                    / self.problem.norbitals**4
                    / forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                    * self.samplers['Worm_P2pp_tau']
                    .postprocess(self.accumulators['Worm_P2pp_tau'].mean)[0, ...]
                ) / mean_sign

            if self.config["QMC"]["WormMeasQQ"]:
                # FIXME: normalization (check/fix)
                result["qqiw-worm"] = (
                    forlib.get_sector_sample_count(get_sector_index(self.config["QMC"]))
                    / forlib.get_sector_sample_count(SECTOR_Z)
                    / self.problem.norbitals**4
                    / forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                    * self.samplers['Worm_QQ_mat']
                    .postprocess(self.accumulators['Worm_QQ_mat'].mean)[0, ...]
                ) / mean_sign

            if self.config["QMC"]["WormMeasQQtau"]:
                # FIXME: normalization (check/fix)
                result["qqtau-worm"] = (
                    forlib.get_sector_sample_count(get_sector_index(self.config["QMC"]))
                    / forlib.get_sector_sample_count(SECTOR_Z)
                    / self.problem.norbitals**4
                    / forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                    * self.samplers['Worm_QQ_tau']
                    .postprocess(self.accumulators['Worm_QQ_tau'].mean)[0, ...]
                ) / mean_sign

            if self.config["QMC"]["WormMeasRaman"]:
                # FIXME: normalization (check/fix)
                result["raman-worm"] = (
                    forlib.get_sector_sample_count(get_sector_index(self.config["QMC"]))
                    / forlib.get_sector_sample_count(SECTOR_Z)
                    / self.problem.norbitals**4
                    / forlib.get_wormeta(get_sector_index(self.config["QMC"]))
                    * self.samplers['Worm_Raman_mat']
                    .postprocess(self.accumulators['Worm_Raman_mat'].mean)[0, ...]
                ) / mean_sign

            if self.worm_config is not None:
                for i, qname in enumerate(self.worm_config['quant_mat']):
                    if qname is None:
                        continue
                    key = 'Worm_Custom_mat{}'.format(i)
                    result[qname] = (
                        forlib.get_sector_sample_count(SECTOR_CUSTOM)
                        / forlib.get_sector_sample_count(SECTOR_Z)
                        / forlib.get_wormeta(SECTOR_CUSTOM)
                        * self.samplers[key]
                        .postprocess(self.accumulators[key].mean)[0, ...]
                    ) / mean_sign

                for i, qname in enumerate(self.worm_config['quant_tau']):
                    if qname is None:
                        continue
                    key = 'Worm_Custom_tau{}'.format(i)
                    result[qname] = (
                        forlib.get_sector_sample_count(SECTOR_CUSTOM)
                        / forlib.get_sector_sample_count(SECTOR_Z)
                        / forlib.get_wormeta(SECTOR_CUSTOM)
                        * self.samplers[key]
                        .postprocess(self.accumulators[key].mean)[0, ...]
                    ) / mean_sign

            self.worm_config = None

        else:  # not worm
            result["sign"] = mean_sign

            result["densitymatrix"] = np.reshape(
                self.accumulators['Z_dm'].mean,
                (2**self.problem.nflavours,
                 2**self.problem.nflavours)
            ) / mean_sign

            result["hist"] = np.reshape(
                self.accumulators['expansion_order'].mean,
                (self.problem.norbitals,
                 2,
                 (self.config["QMC"]["MaxHisto"] + 1))
            )

            result["ntau-n0"] = np.reshape(
                self.accumulators['Z_ntau_n0'].mean,
                (self.problem.norbitals,
                 2,
                 self.problem.norbitals,
                 2,
                 (self.config["QMC"]["Ntau"] + 1))
            ) / mean_sign

            result["gtau-full"] = np.reshape(
                self.samplers['Z_G_tau']
                .postprocess(self.accumulators['Z_G_tau'].mean),
                (self.problem.norbitals,
                 2,
                 self.problem.norbitals,
                 2,
                 -1)
            ) / mean_sign

            result["gleg-full"] = np.reshape(
                self.samplers['Z_G_leg']
                .postprocess(self.accumulators['Z_G_leg'].mean),
                (self.problem.norbitals,
                 2,
                 self.problem.norbitals,
                 2,
                 -1)
            ) / mean_sign

            # XXX figure out why we have a sign here
            result["giw-full"] = -1.0 * (
                np.reshape(
                    self.samplers['Z_G_mat']
                    .postprocess(self.accumulators['Z_G_mat'].mean),
                    (self.problem.norbitals,
                     2,
                     self.problem.norbitals,
                     2,
                     -1))
            ) / mean_sign
            
            if self.config["QMC"]["MeasG2iw"]:
                result["g2iw-full"] = -1.0 * (
                    np.reshape(
                        self.samplers['Z_G2iw_mat']
                        .postprocess(self.accumulators['Z_G2iw_mat'].mean),
                        (self.problem.norbitals,
                        2,
                        self.problem.norbitals,
                        2,
                        2*self.config["QMC"]["Niw"],
                        2*self.config["QMC"]["Niw"]))
                ) / mean_sign

            if self.config["QMC"]["FourPnt"]==8:
                result['g4iw-full'] = np.reshape(
                    self.accumulators['Z_G4iw_full'].mean,
                    (self.problem.norbitals,
                    2,
                    self.problem.norbitals,
                    2,
                    self.problem.norbitals,
                    2,
                    self.problem.norbitals,
                    2,
                    -1), order="F"
                ) / mean_sign

            ### TODOcompl: for some unknown reason these
            # quantities have to be complex conjugated
            result["gleg-full"] = np.conj(result["gleg-full"])
            result["gtau-full"] = np.conj(result["gtau-full"])

        return result

    def get_resultset(self, qmc_config):
        # Build up result array from QMC data
        # FIXME TODO: all quantities of interest that could be
        # obtained in the old version should be ported to the new
        # version
        result = {
            # "rhist": ctqmc.rhisto,
            # "lhist": ctqmc.lhisto,
            "ssts-states": forlib.get_FIXME("StatesSuperstates",
                                            (2**self.problem.nflavours),
                                            np.int32),
            "energies-eigenstates": forlib.get_FIXME("EigenstatesEnergies",
                                                     (2**self.problem.nflavours),
                                                     np.double),
            "occbasis-mapping": forlib.get_FIXME("OccbasisMapping",
                                                 (2**self.problem.nflavours,
                                                  self.problem.norbitals, 2),
                                                 np.int32),
            # "hist-sst": ctqmc.outersuperstatehisto,
            # "hist-state": ctqmc.outerstatehisto,
            # "sign-sst": ctqmc.signhistosuperstates,
            # "sign-state": ctqmc.signhistostates,
            # "contrib-sst": ctqmc.tracecontribsuperstates,
            # "contrib-state": ctqmc.tracecontribstates,
            "occ": forlib.get_FIXME("occ",
                                    (self.problem.norbitals, 2,
                                     self.problem.norbitals, 2),
                                    np.double),
            "accept-ins": forlib.get_FIXME("AccAdd",
                                           (MAX_SECTOR),
                                           np.double),
            "accept-ins4": forlib.get_FIXME("AccAdd4",
                                            1,
                                            np.double)[0],
            "accept-rem": forlib.get_FIXME("AccRem",
                                           (MAX_SECTOR),
                                           np.double),
            "accept-rem4": forlib.get_FIXME("AccRem4",
                                            1,
                                            np.double)[0],
            # "accept-pair-tau": ctqmc.accpair,
            "accept-glob": forlib.get_FIXME("AccGlob",
                                            1,
                                            np.double)[0],
            "accept-shift": forlib.get_FIXME("AccShift",
                                             1,
                                             np.double)[0],
            "accept-flavourchange": forlib.get_FIXME("AccFlavc",
                                                     1,
                                                     np.double)[0],
            "time-warmup": forlib.get_FIXME("time_warmup",
                                            1,
                                            np.double)[0],
            "time-simulation": forlib.get_FIXME("time_sim",
                                                1,
                                                np.double)[0],
            "time-sampling": forlib.get_FIXME("time_sim_steps",
                                              1,
                                              np.double)[0],
        }
        # if qmc_config["MeasSsquare"] != 0:
        #     result["ssquare-operator"] = ctqmc.ssquare_operator
        # if qmc_config["Gtau_mean_step"] != 0:
        #     result["gtau-mean-step"] = ctqmc.gtau_mean_step
        # if qmc_config["Gtau_mid_step"] != 0:
        #     result["gtau-mid-step"] = ctqmc.gtau_mid_step
        # if qmc_config["sign_step"] != 0:
        #     result["sign-step"] = ctqmc.sign_step
        # if qmc_config["MeasGSigmaiw"] != 0:  # enable if GSigma Z-meas fixed
        #     result["gsigmaiw"] = ctqmc.gsigmaiw
        # if qmc_config["MeasG2iw"] != 0:
        #     result["g2iw"] = ctqmc.g2iw
        if qmc_config["MeasDensityMatrix"] != 0:
            result["rho2"] = forlib.get_FIXME(
                "rho2",
                (self.problem.nflavours,) * 4,
                np.double
            ).reshape((self.problem.norbitals, 2) * 4)
            result["rho1"] = forlib.get_FIXME(
                "rho1",
                (self.problem.nflavours,) * 2,
                np.cdouble
            ).reshape((self.problem.norbitals, 2) * 2)
        # if qmc_config["MeasSsquare"] != 0 and (qmc_config["MeasDensityMatrix"] != 0
        #                                        or qmc_config["Eigenbasis"] != 0):
        #     result["s2"] = ctqmc.ssquare
        # if qmc_config["MeasExpResDensityMatrix"] != 0:
        #     result["expresdensitymatrix"] = ctqmc.expresdensitymatrix
        # if qmc_config["FourPnt"] & 1:
        #     result["g4tau"] = ctqmc.g4tau
        # if qmc_config["FourPnt"] & 2:
        #     result["g4leg"] = ctqmc.g4leg
        # if qmc_config["FourPnt"] & 4:
        #     result["g4iw"] = ctqmc.g4iw
        # if qmc_config["FourPnt"] & 4 and qmc_config["MeasG4iwPP"] != 0:
        #     result["g4iw-pp"] = ctqmc.g4iw_pp
        # if qmc_config["segment"] != 0:
        #     result["hist-seg"]= ctqmc.histo_seg

        return result


    @classmethod
    def get_resultset_worm(cls, qmc_config):
        # Build up result array from QMC data
        result = {
            "accept-worm-ins": forlib.get_FIXME("AccWormAdd",
                                                (MAX_SECTOR-1),
                                                np.double),
            "accept-worm-rem": forlib.get_FIXME("AccWormRem",
                                                (MAX_SECTOR-1),
                                                np.double),
            "accept-worm-rep": forlib.get_FIXME("AccWormRep",
                                                (MAX_SECTOR-1),
                                                np.double),
            "steps-worm-partition": np.array([forlib.get_sector_sample_count(sector)
                                              for sector in range(1, MAX_SECTOR + 1)],
                                             np.double),
            "worm-eta": np.array([forlib.get_wormeta(sector)
                                  for sector in range(1, MAX_SECTOR + 1)]),
            # The conversion to double for large integers is necessary because
            # otherwise squaring it for the calculation of errors produces
            # negative values and nan errorbars.
        }

        return result
