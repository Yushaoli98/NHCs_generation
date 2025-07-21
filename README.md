# NHC_AuCl_Generator

用于从NHCs分子的SMILES结构中，自动识别卡宾碳，并拼接 AuCl 基团形成配合物，再使用 XTB 对三维结构进行优化。



| `config/settings.py` | 全局配置参数（路径、键长、XTB路径等） |
| `nhc_generator/detect.py` | 识别碳负离子，即卡宾中心 |
| `nhc_generator/conformer.py` | UFF 多构象搜索并选出最低能量构象 |
| `nhc_generator/geometry.py` | 分子三维坐标提取、方向选择、Van der Waals 半径获取等 |
| `nhc_generator/aucl_attach.py` | 生成Au和Cl原子坐标 |
| `nhc_generator/io_utils.py` | 读写XYZ结构文件 |
| `nhc_generator/xtb_wrapper.py` | 调用XTB进行结构优化（支持并行） |
| `nhc_generator/process.py` | 封装整个SMILES → NHC → NHC-AuCl → 优化的流水线 |
| `main.py` | 程序主入口，支持多进程执行 |

## 📦 安装依赖

```bash
conda create -n nhc_gen python=3.9
conda activate nhc_gen
conda install -c conda-forge rdkit scipy pandas numpy
