# 1D2V Fluid Equation Solver

这是一个一维空间、二维速度分量的电磁场-电子流体耦合求解器。程序模拟左边界入射激光与有限长度等离子体片相互作用的过程，沿 `z` 方向离散空间，同时演化横向速度 `v_x` 和纵向速度 `v_z`，因此可称为 1D2V Maxwell-Fluid 模型。

## 求解的方程组

程序使用无量纲变量求解一维 Maxwell 方程与电子流体方程的耦合系统。空间坐标为 `z`，电磁场保留 `E_x`、`E_z` 和 `B_y`，电子流体变量为电子数密度 `n_e`、横向速度 `v_x` 和纵向速度 `v_z`：
```
$$
\begin{cases}
    -\dfrac{\partial \tilde{B}_{y}}{\partial \tilde{z}}
    =
    -\tilde{n}_{e}\tilde{v}_{x}
    +
    \dfrac{\partial \tilde{E}_{x}}{\partial \tilde{t}},
    & \text{Ampere 定律（横向电场）} \\[8pt]

    \dfrac{\partial \tilde{E}_{x}}{\partial \tilde{z}}
    =
    -\dfrac{\partial \tilde{B}_{y}}{\partial \tilde{t}},
    & \text{Faraday 定律（横向磁场）} \\[8pt]

    \tilde{n}_{e}\tilde{v}_{z}
    =
    \dfrac{\partial \tilde{E}_{z}}{\partial \tilde{t}},
    & \text{纵向总电流平衡} \\[8pt]

    \dfrac{\partial \tilde{n}_{e}}{\partial \tilde{t}}
    +
    \dfrac{\partial(\tilde{n}_{e}\tilde{v}_{z})}{\partial \tilde{z}}
    =
    0,
    & \text{连续性方程} \\[8pt]

    \dfrac{\partial \tilde{v}_{x}}{\partial \tilde{t}}
    +
    \tilde{v}_{z}\dfrac{\partial \tilde{v}_{x}}{\partial \tilde{z}}
    =
    -\left(\tilde{E}_{x}-\tilde{v}_{z}\tilde{B}_{y}\right),
    & \text{横向动量方程} \\[8pt]

    \dfrac{\partial \tilde{v}_{z}}{\partial \tilde{t}}
    +
    \tilde{v}_{z}\dfrac{\partial \tilde{v}_{z}}{\partial \tilde{z}}
    =
    -\left(\tilde{E}_{z}+\tilde{v}_{x}\tilde{B}_{y}\right)
    -
    \beta \dfrac{1}{\tilde{n}_{e}}
    \dfrac{\partial \tilde{n}_{e}}{\partial \tilde{z}},
    & \text{纵向动量方程}
\end{cases}
$$

其中压力项系数

$$
\beta = \frac{\gamma_e k_B T_e}{m_e c^2}.
$$
```
归一化关系为：

- 时间：`\tilde{t} = \omega_{pe} t`
- 空间：`\tilde{z} = \omega_{pe} z / c`
- 速度：`\tilde{v} = v / c`
- 密度：`\tilde{n}_e = n_e / n_{e0}`
- 电场：`\tilde{E} = E / (m_e c \omega_{pe} / e)`
- 磁场：`\tilde{B} = B / (m_e \omega_{pe} / e)`
- 电子等离子体频率：`\omega_{pe} = \sqrt{n_{e0} e^2 / (\epsilon_0 m_e)}`

代码中 `qz` 当前等同于纵向无量纲速度 `v_z`，`j_x = n_e v_x`，`j_z = n_e v_z`。

## 数值模型

- 空间维度：一维 `z` 方向。
- 速度维度：电子流体速度包含 `v_x` 与 `v_z` 两个分量。
- 电磁场：演化 `E_x`、`E_z`、`B_y`。
- 激光边界：默认从左边界注入，电场为正弦驱动，可选择 `sin^2` 平滑上升包络。
- 等离子体区域：默认平台区间为 `[5 um, 20 um]`，两侧可设置密度缓变层。
- 边界与界面：支持真空-等离子体界面跳跃条件开关；默认关闭。
- 特殊测试模式：`linear_overdense_test` 可冻结 `n`、`E_z`、`v_z` 和 `q_z`，只推进 `E_x`、`B_y` 和 `v_x`。

## 项目结构

```text
.
├── CMakeLists.txt
├── README.md
├── src/
│   ├── main.cpp
│   ├── config.h
│   ├── fluid_maxwell_1d.cpp
│   ├── fluid_solver.cpp
│   ├── maxwell_solver.cpp
│   ├── result_output.cpp
│   └── physics_constants.h
├── PostProcessing/
│   ├── plot_output.py
│   ├── plot_time_evolution_at_positions.py
│   ├── plot_spatial_profile_for_state.py
│   └── config.py
├── output/
└── results/
```

## 编译

项目使用 CMake 构建，要求支持 C++11 的编译器。

```bash
cmake -S . -B build
cmake --build build
```

Windows 下生成的可执行文件通常位于：

```text
build/fluid_maxwell_1d.exe
```

Linux 或 macOS 下通常为：

```text
build/fluid_maxwell_1d
```

## 运行

```bash
./build/fluid_maxwell_1d
```

程序启动后会打印网格、时间步长、等离子体频率、激光参数、等离子体区间等信息，并将结果写入 `output/`。

## 主要参数

主要物理和数值参数在 `src/config.h` 的 `Config` 结构体中设置：

| 参数 | 含义 | 默认值 |
| --- | --- | --- |
| `lambda0` | 入射激光波长 | `3.0e-6 m` |
| `intensity_w_cm2` | 激光强度 | `1.0e17 W/cm^2` |
| `electron_temperature_ev` | 电子温度 | `100 eV` |
| `electron_density0` | 参考电子密度 | `3.0e26 m^-3` |
| `length` | 计算区域长度 | `25.0e-6 m` |
| `dz` | 空间步长 | `0.002e-6 m` |
| `t_end` | 模拟终止时间 | `120 fs` |
| `snapshot_dt` | 输出快照时间间隔 | `0.5 fs` |
| `plasma_left` / `plasma_right` | 等离子体平台区间 | `5 um` / `20 um` |
| `plasma_ramp_width` | 等离子体边缘缓变层宽度 | `0.5 um` |
| `laser_from_left_boundary` | 是否从左边界注入激光 | `true` |
| `laser_smooth_ramp` | 是否使用平滑激光包络 | `false` |

## 输出文件

运行结果默认写入 `output/`：

- `state_0000000.csv`、`state_*.csv`：不同时间步的空间剖面。
- `diagnostics.csv`：整体诊断量时间序列。

`state_*.csv` 的列包括：

| 列名 | 含义 |
| --- | --- |
| `z` | 空间位置，单位 m |
| `Ex` | 横向电场，单位 V/m |
| `Ez` | 纵向电场，单位 V/m |
| `By` | 横向磁场，单位 T |
| `ne` | 电子数密度，单位 m^-3 |
| `vx` | 横向电子速度，单位 m/s |
| `vz` | 纵向电子速度，单位 m/s |
| `Pe` | 电子压力，单位 Pa |

`diagnostics.csv` 的列包括 `time_fs`、`n_min`、`n_max`、`vx_rms`、`vz_rms`、`energy_em` 和 `energy_fluid`。

## 后处理

`PostProcessing/` 中提供了绘图脚本，可绘制指定位置处的时间演化或某个快照的空间剖面。

```bash
python PostProcessing/plot_output.py time --positions-um 0.5,2.0,4.0
python PostProcessing/plot_output.py space --state-file latest
```

如需使用归一化输出，可增加 `--normalize`：

```bash
python PostProcessing/plot_output.py space --state-file latest --normalize
```

后处理脚本默认从 `output/` 读取数据，并将图像写入 `results/`。
