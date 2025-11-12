cmake_minimum_required(VERSION 3.9...3.31)

function(declare_compileall_test test_name src_dir)
    add_test(NAME ${test_name}
        COMMAND ${Python3_EXECUTABLE} -m compileall -q ${src_dir})
endfunction()

function(declare_regression_test test_name test_dir test_script)
    add_test(
        NAME ${test_name}
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/regressiontest.py
            ${PROJECT_SOURCE_DIR} ${test_dir} ${test_script}
        )
    set_tests_properties(${test_name}
	    PROPERTIES ENVIRONMENT PYTHON_EXECUTABLE=${Python3_EXECUTABLE}
        )
endfunction()

function(add_test_aux_file filename)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
        ${CMAKE_CURRENT_BINARY_DIR}/${filename}
        COPYONLY)
endfunction()


enable_testing()
declare_compileall_test(1-syntax_check ${PROJECT_SOURCE_DIR})

add_test_aux_file(bsp_444/Parameters.in)
add_test_aux_file(bsp_444/e0)
add_test_aux_file(bsp_444/ftau.dat)
add_test_aux_file(bsp_444/fiw.dat)
add_test_aux_file(bsp_444/reference_data.hdf5)
add_test_aux_file(bsp_444/axes.dat)
add_test_aux_file(bsp_444/quantities.dat)

declare_regression_test(2-regressiontest-bsp_444 bsp_444 cthyb)

add_test_aux_file(bsp_444___compl_vk___hopping_compl/Parameters.in)
add_test_aux_file(bsp_444___compl_vk___hopping_compl/e0)
add_test_aux_file(bsp_444___compl_vk___hopping_compl/ftau.dat)
add_test_aux_file(bsp_444___compl_vk___hopping_compl/fiw.dat)
add_test_aux_file(bsp_444___compl_vk___hopping_compl/reference_data.hdf5)
add_test_aux_file(bsp_444___compl_vk___hopping_compl/axes.dat)
add_test_aux_file(bsp_444___compl_vk___hopping_compl/quantities.dat)

declare_regression_test(3-regressiontest-bsp_444_complex bsp_444___compl_vk___hopping_compl cthyb)

add_test_aux_file(strontium_iridate/Parameters.in)
add_test_aux_file(strontium_iridate/wannier90_hk_t2gbasis.dat_)
add_test_aux_file(strontium_iridate/one_more_iter-2022-10-11-Tue-14-00-53.hdf5)
add_test_aux_file(strontium_iridate/reference_data.hdf5)
add_test_aux_file(strontium_iridate/axes.dat)
add_test_aux_file(strontium_iridate/quantities.dat)

add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0/axes.dat)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0/epsk)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0/Parameters.in)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0/quantities.dat)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0/reference_data.hdf5)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0/Vk)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0/w2dlevels.dat)

declare_regression_test(4-regressiontest-2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0 2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0 cthyb)

add_test_aux_file(3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00/axes.dat)
add_test_aux_file(3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00/epsk)
add_test_aux_file(3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00/Parameters.in)
add_test_aux_file(3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00/quantities.dat)
add_test_aux_file(3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00/reference_data.hdf5)
add_test_aux_file(3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00/Vk)
add_test_aux_file(3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00/w2dlevels.dat)

declare_regression_test(5-regressiontest-3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00 3orb_bws2.5_3.0_3.5_lvls1.0_0.5_0.0_kanU1.6_J0.4_mu5.00 cthyb)

add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm/axes.dat)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm/epsk)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm/Parameters.in)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm/quantities.dat)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm/reference_data.hdf5)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm/Vk)
add_test_aux_file(2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm/w2dlevels.dat)

declare_regression_test(6-regressiontest-2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm 2orb_bws2.0_4.0_lvls0.0_kanU8.0_J2.0_mu7.0_worm cthyb)

add_test_aux_file(1orb_bethe_dmft/axes.dat)
add_test_aux_file(1orb_bethe_dmft/Parameters.in)
add_test_aux_file(1orb_bethe_dmft/quantities.dat)
add_test_aux_file(1orb_bethe_dmft/reference_data.hdf5)

declare_regression_test(7-regressiontest-1orb_bethe_dmft 1orb_bethe_dmft DMFT.py)
