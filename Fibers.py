# rebarY = np.array([-200, -100, 100, 200])
# rebarZ = np.array([-150, 0, 150])
# rebarYZ = np.zeros([len(rebarY) * len(rebarZ), 2])
# for ii, Y in enumerate(rebarY):
#     for jj, Z in enumerate(rebarZ):
#         rebarYZ[ii*len(rebarZ) + jj, :] = [Y, Z]
# for YZ in rebarYZ:
#     fiber(*YZ, Ay, 1)

# import openseespy.postprocessing.ops_vis as opsv
# uniaxialMaterial('Steel01', 4, 355, 200000, 0.02)
# Mat = 4
# fib_sec_1 = [['section', 'Fiber', Mat, '-GJ', 1.0e6],
#              ['patch', 'quad', Mat, 4, 1,  0.032, 0.317, -0.311, 0.067, -0.266, 0.005, 0.077, 0.254],
#              ['patch', 'quad', Mat, 1, 4,  -0.075, 0.144, -0.114, 0.116, 0.075, -0.144, 0.114, -0.116],
#              ['patch', 'quad', Mat, 4, 1,  0.266, -0.005,  -0.077, -0.254,  -0.032, -0.317,  0.311, -0.067]
#              ]

# opsv.fib_sec_list_to_cmds(fib_sec_1)
# matcolor = ['r', 'lightgrey', 'gold', 'blue', 'w', 'w']
# opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor)
# plt.axis('equal')
# plt.show()