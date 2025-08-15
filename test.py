# # import cmath, math
# # x1 = [1, 0, 0, 0]
# # x2 = [0, 1, 0, 0]
# # N = len(x1)

# # def dft(x):
# #     return [sum(x[n] * cmath.exp(-2j * math.pi * k * n / N)
# #                 for n in range(N)) for k in range(N)]

# # X1 = dft(x1)
# # X2 = dft(x2)

# # print('X1 magnitudes:', [round(abs(v), 4) for v in X1])
# # print('X2 magnitudes:', [round(abs(v), 4) for v in X2])
# # print('Are magnitudes equal?', all(abs(a) - abs(b) < 1e-9
# #                                    for a, b in zip(X1, X2)))

# import cmath, math

# N = 8
# # 构造简单图像：中央 4x4 方块为 1
# img = [[0]*N for _ in range(N)]
# for i in range(2, 6):
#     for j in range(2, 6):
#         img[i][j] = 1

# def dft2d(x):
#     N = len(x)
#     return [[sum(x[m][n] * cmath.exp(-2j*math.pi*(k*m + l*n)/N)
#                  for m in range(N) for n in range(N))
#              for l in range(N)] for k in range(N)]

# def idft2d(X):
#     N = len(X)
#     return [[sum(X[k][l] * cmath.exp(2j*math.pi*(k*m + l*n)/N)
#                  for k in range(N) for l in range(N)) / (N*N)
#              for n in range(N)] for m in range(N)]

# Y = dft2d(img)                            # 完整 k-space
# mask = [[1 if k % 2 == 0 else 0 for l in range(N)]
#         for k in range(N)]               # 每隔一行采样
# Y_us = [[Y[k][l]*mask[k][l] for l in range(N)] for k in range(N)]
# recon = idft2d(Y_us)                     # 零填充逆变换

# # 打印原图与欠采样重建的 8x8 矩阵（四舍五入）
# print("Original image:")
# for row in img:
#     print(row)

# print("\nReconstructed from undersampled k-space:")
# for row in recon:
#     print([round(r.real, 2) for r in row])

import h5py
import numpy as np
import matplotlib.pyplot as plt

# 读取一个样本
with h5py.File('D:/MyDatasets/fastMRI/knee_singlecoil_test~/singlecoil_test/file1000174.h5', 'r') as f:
    kspace = f['kspace'][()]          # 复数 k-space, shape: (num_slices, height, width)
    print('k-space shape:', kspace.shape, 'dtype:', kspace.dtype)

    # 查看数值范围
    print('k-space min/max:', kspace.min(), kspace.max())

    # 取第0张切片做逆FFT，得到图像
    img = np.fft.ifft2(kspace[30], norm='ortho')
    magnitude = np.abs(img)

plt.imshow(magnitude, cmap='gray')
plt.title('Slice 0 magnitude image')
plt.axis('off')
plt.show()


import numpy as np
from py3dgespar import run_gespar3d_stats

# 假设 kspace_stacked.shape = (num_slices, height, width)
# 先做 2D IFFT 得到复图像，然后堆叠成 3D 体 volume
volume = np.fft.ifft2(kspace_stacked, norm='ortho')
volume = np.abs(volume)               # 仅保留幅度，模拟“缺相位”场景

# 计算其 3D FFT 幅度平方作为测量 y
y = np.abs(np.fft.fftn(volume))**2

# 设置稀疏度 k、维度 dimlen 等
dimlen = volume.shape[0]              # 假设立方体，需自行裁剪或插值
k = 20                                # 稀疏度
m = y.size                            # 测量数量
max_t = 1000
snr = np.inf

nmse, supp, t_arr, iter_arr = run_gespar3d_stats(
    volume.flatten(), dimlen, k, m, max_t, snr, verbose=True
)
print('NMSE:', nmse)
