import numpy as np
import w2dyn.dmft.interaction as interaction

norb = 3

i = interaction.Kanamori(norb, 2.0, 0.4, 0.8)
np.savetxt("u_matrix.dat", np.ravel(i.u_matrix))

m = 0.2 * np.eye(norb)
m[1, 1] += 0.03
m[2, 2] += 0.05
m = m[:, np.newaxis, :, np.newaxis] * np.eye(2)[np.newaxis, :, np.newaxis, :]
np.savetxt("muimp.dat", np.ravel(m))
