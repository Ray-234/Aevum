<p align="center">
  <img src="docs/assets/aevum-logo.svg" alt="Aevum logo" width="760">
</p>

# aevum — 行星深时演化引擎 + 世界历史档案 + 策略地图编译器

[English README](README_EN.md) · [文档索引](docs/INDEX.md) · [结果展示](docs/RESULT_SHOWCASE.md)

`aevum`（拉丁语“纪元/世”）不是一个“地图生成器”，而是一个能回答
**“这里为什么会变成这样”** 的因果系统。它从行星形成一路演化到“现代”，
保存完整的形成历史，并把真值世界**编译**成一张《文明》式的六边格策略地图。

> 核心理念（与规划一致）：**先把每个特性的“契约”预留好，再逐个填补占位模块。**
> 预留的不是几百个空字段，而是每个特性的*定义、单位、依赖、时空表示、生产模块、
> 保真度等级、不确定性与验证方法*。空模块也能串起来跑通整条因果链。

当前为 **v0.1**：全链路因果完整、能在数秒内跑完整个深时历史并产出现代策略地图与
因果时间线。每替换一个占位模块，世界会更可信，但架构无需推倒重来。

## 当前交接入口（2026-07-08）

后续在 CPU 集群继续实验时，先读
[`docs/PROJECT_HANDOFF_20260708.md`](docs/PROJECT_HANDOFF_20260708.md) 与
[`docs/CLUSTER_EXPERIMENT_PLAN_20260708.md`](docs/CLUSTER_EXPERIMENT_PLAN_20260708.md)。
当前决策是：板块/地形生成阶段暂时收官；内置快速气候引擎冻结为诊断原型；
下一阶段先用真实地球子图与外部气候工具校准，再把月平均气温/降水后处理成
Köppen 与 biome。

## 结果展示

更多图和来源说明见
[`docs/RESULT_SHOWCASE.md`](docs/RESULT_SHOWCASE.md)。README 中只保留少量真实输出，
不提交本地完整 `out*` 实验目录。

![72000-cell elevation snapshot](docs/assets/showcase/elevation_72000_seed707.png)

### 最新地形演化视频

<video src="docs/assets/showcase/elevation_evolution_earthlike_seed42.mp4" controls width="100%" poster="docs/assets/showcase/elevation_72000_seed707.png"></video>

如果当前 Markdown 渲染器不内嵌播放视频，可以直接打开
[`elevation_evolution_earthlike_seed42.mp4`](docs/assets/showcase/elevation_evolution_earthlike_seed42.mp4)。
该视频来自 `out_elevation_evolution_videos_6worlds_20260706/earthlike_seed42/`，
展示类地世界从深时构造演化到终端地形图的高程变化。

| 构造/地貌诊断 | 洋底与造山语义 |
|---|---|
| ![Terrain diagnostic contact sheet](docs/assets/showcase/terrain_contact_sheet_72000_seed707.png) | ![Orogenic hierarchy overlay](docs/assets/showcase/orogenic_hierarchy_72000_seed707.png) |

| 海底地貌 | 构造对象层 |
|---|---|
| ![Bathymetry classes](docs/assets/showcase/bathymetry_72000_seed707.png) | ![Tectonic object layer](docs/assets/showcase/tectonic_objects_earthlike_seed42.png) |

> 气候图暂不作为最终产品展示。仓库中保留的 temperature / precip / biome 图只代表旧的
> 内置快速气候原型；当前研发路线是先用真实地球子图与外部气候工具校准，再生成
> Köppen 与 biome。

---

## 三层严格分离

| 层 | 模块 | 说明 |
|---|---|---|
| **真值层** | `WorldState` | 尽量物理一致的行星状态（场/网络/对象/全局标量） |
| **历史层** | `WorldArchive` + `EventBus` | 事件、谱系、地层、矿床、环境变化与因果链 |
| **游戏层** | `MapCompiler` | 压缩成海洋/平原/丘陵/山脉/河流/资源/地块产出 |

游戏平衡**只**发生在第三层，绝不反向篡改真值层。内部时间轴是“自行星形成以来的
时间（Myr）”，“冥古宙/第四纪”等只是地球预设场景的显示标签。

## 架构总览

```
PlanetSpec → FeatureRegistry → WorldState → DeepTimeScheduler
   → 全历史占位流程（stellar→interior→impacts→tectonics→terrain
                      →climate→biogeochem→biosphere→resources）
   → WorldArchive → MapCompiler → 逐模块替换与验证
```

## 快速开始

```bash
python3.13 -m venv .venv && ./.venv/bin/pip install -e .
# 跑一颗类地行星，输出地图/图层/时间线到 out/
./.venv/bin/python -m aevum.cli run --preset earthlike --cells 8000 --out out
# 列出 6 个基准世界
./.venv/bin/python -m aevum.cli presets
# 检查特性注册表（契约 + 反馈环）
./.venv/bin/python -m aevum.cli registry
```

输出包含：`elevation/temperature/precip/biomes/plates/crust_age` 图层、
`hexmap.png`（编译后的六边格地图：河流、资源、起点）、`timeline.png`
（深时事件时间线）、`history.png`（古地理演化）、以及
`timeline.json / explain_examples.json / lineages.json / spec.json`。

## 六个基准世界（`aevum/spec/presets.py`）

`earthlike`(移动盖层) · `waterworld`(水世界) · `arid`(干旱) ·
`stagnant_lid`(停滞盖层) · `tidally_locked`(潮汐锁定) · `frozen`(冰封)。
差异来自三个层面：**参数差异**（恒星/轨道/组成/大气）、**体制差异**
（`TectonicRegime`：移动/缓慢/间歇/停滞盖层）、**历史差异**（撞击/火山/裂谷/
创新/灭绝产生的分叉历史）。Earth 只是其中一个测试场景，不是唯一校准目标。

## 规划 → 代码映射

| 规划要素 | 实现位置 |
|---|---|
| Feature Registry（特性契约） | `aevum/core/registry.py` + `aevum/features.py`（53 条） |
| 四类世界数据（场/网络/对象/事件） | `aevum/core/state.py`、`events.py` |
| 统一模块接口 `step(...)→delta,events,diag` | `aevum/core/module.py` |
| 可重复随机数（种子+模块+时间+事件） | `aevum/core/rng.py` |
| 球面分层网格（DGGS 思路） | `aevum/core/grid.py`（Fibonacci+球面Voronoi） |
| 多速率/事件驱动/准平衡调度 | `aevum/core/scheduler.py`（气候按阈值再求解） |
| 恒星与轨道 / 内部热史 | `modules/stellar.py` / `modules/interior.py` |
| 构造地质（球面板块+事件规则） | `modules/tectonics.py`（拉格朗日平流栅格化） |
| 地外撞击（时变随机+尺度律） | `modules/impacts.py` |
| 地形/侵蚀/沉积/海平面/河网 | `modules/terrain.py` + `modules/hydrology.py` |
| 大气海洋（EBM+风带+雨影+海冰） | `modules/climate.py` |
| 生地化（COPSE 式 C/O 箱循环+风化反馈） | `modules/biogeochem.py` |
| 生物圈（功能类群/生态位/扩散/创新/灭绝） | `modules/biosphere.py` |
| 资源成因谱系（形成→保存→埋深→品位） | `modules/resources.py` |
| 策略地图编译（重采样→分类→河流→港湾→资源→产出→起点公平） | `compiler/map_compiler.py` |
| 五类验证（守恒/拓扑/因果/回归/世界组合） | `aevum/validation.py` + `tests/` |

## 验证体系（“不能只看起来像”）

- **守恒**：碳库存（大气↔深部）精确闭合（相对误差 ~1e-16）。
- **拓扑**：板块铺满球面无空洞、河网无环、入海连通。
- **因果**：每个矿床有成因/年代/围岩；每座高地可溯源（造山/火山/厚地壳）。
- **回归**：固定种子逐位复现。
- **世界组合**：6 个基准世界均通过守恒/拓扑测试。

每个输出附带 provenance（数值/单位/生产模块/保真度/不确定区间/直接原因/上游事件），
点任意一格即可回答“为什么是山 / 为什么这里有铜没有石油 / 这条河为何这样流”。

## “地球史机制族”——预留的是机制，不是地球专名

引擎不硬编码“寒武纪大爆发”，而预留能在别的星球上产生类似/不同/完全不发生
之结果的通用机制：构造体制转换、超大陆循环（板块周期性重组）、海洋通道、
长期碳—硅酸盐循环（气候恒温器）、氧化还原跃迁（大氧化事件自发涌现）、
冰—反照率反馈、营养盐—生产力—缺氧、生态创新（光合/多细胞/登陆/森林）、
灾变与恢复（撞击/灭绝），以及资源成因谱系与（未来的）技术圈反馈。

## 路线图（阶段 0→6）

- **阶段 0（已完成）**：契约 + 骨架（规范/注册表/网格/接口/调度/档案/RNG/6 世界）。
- **阶段 1（已完成）**：低保真全历史薄切片，产出现代六边格地图 + 事件时间线。
- **阶段 2+（占位/待提升）**：构造可信度（应力/流变/真实板块拓扑、可借鉴 GPlates 数据模型）、
  气候海洋（季节/三维环流，可用 ExoPlaSim 关键帧校验）、生地化（向 cGENIE 靠拢）、
  地表过程（借鉴 Badlands 的过程划分）、资源概率化品位—吨位、地图编辑器与全球三维浏览。
  注册表中仍为 `reserved`/低保真的特性即后续填补目标。

## 设计取舍

不复制地球事件，复制产生事件的机制；不只存地图状态，存其形成历史；
不用噪声直接生成山脉/沙漠/矿产（噪声只用于初始扰动与随机事件）；
模块不共享单一时间步与空间网格；机器学习只作代理/反演/筛选，不作因果核心；
渲染放在最后。
