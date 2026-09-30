# 环状 DNA 条码拼接服务（barcode-assembler）

一批来自同一环状条码的等长 DNA 短序列，判定能否拼成**可信的唯一条码**，
避免把重复片段在不同排列下的差异误报为不同病原株。

## 方法

每条长度 `L = k+1` 的读数视为德布鲁因（de Bruijn）多重图中的一条有向边
`s[:-1] → s[1:]`；**每次出现都是一条独立平行边（一份证据）**。
"每条读数恰好使用一次、相邻重叠长度为 `L-1` 的闭环"等价于该多重图的
欧拉回路。

- **唯一**：所有闭环拼法同属一个规范等价类；
- **歧义**：存在两个以上规范等价类，返回两条不同规范见证；
- **无解**：给出可复核原因——`degree_imbalance`（入度≠出度的 k-mer 清单）
  和/或 `fragmented_graph`（非零度顶点的弱连通分量）。

规范等价归一化：

1. 环的循环移位；
2. 整条环的反向互补；
3. 相同读数（平行边）互换不产生新类别（枚举时同标签只取当前最小输入序号，
   因而**使用次序按输入序号稳定分配**且跨次运行确定）。

枚举从固定起点出发做带桥剪枝的欧拉回路 DFS，按规范代表去重，找到 2 个
不同等价类即停。约束（6~24 条、长度 3~8）下毫秒级完成。

**方向未知模式**（请求加 `"unknown_orientation": true`）下，每条读数可独立
取原序列或反向互补：在"每个相同标签取多少份反向互补"的变量空间上联合回溯
定向与闭环次序（增量维护各 k-mer 出度-入度差，用剩余变量的逐点贡献上下界
剪枝），再做连通检查与欧拉枚举；整批翻转的对称性通过固定一个非回文锚定读数
的定向消除，结果稳定。详见下文 API 小节。

## API

零第三方依赖（Python 标准库）。

- `GET /health` → `{"status":"ok",...}`
- `POST /assemble`，请求体：

```json
{ "sequences": ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"] }
```

响应 `status` 为 `unique` / `ambiguous` / `no_solution`（非法输入返回 400）。

`unique` 响应包含：

- `canonical_barcode`：规范条码（环旋转+反向互补下的字典序最小代表）；
- `barcode`：与 `order` 同方向的环；
- `order`：使用次序，元素为 **1 起始的输入序号**（每个序号恰好出现一次）；
- `evidence`：相邻重叠证据（`prev`/`next` 序号、`overlap`、`overlap_length = L-1`、
  `appended_base`），最后一条与首条首尾相接。

`ambiguous` 响应在 `witnesses` 中给出两条不同规范见证（结构同上）。

### 方向未知模式（`unknown_orientation`）

部分测序批次会丢失单条读数的链方向。在既有请求中加入布尔标志即可启用，
**不启用时输入、结论与证据结构与旧版完全一致**：

```json
{
  "sequences": ["ATT", "ATC", "TCG", "GCG", "TGC", "TTG"],
  "unknown_orientation": true
}
```

启用后：

- 每条读数【独立】可按原序列（`forward`）或其反向互补
  （`reverse_complement`）参与拼接，但每个输入序号仍然**恰好使用一次**；
- **所有逐条定向与闭环次序联合求解**（在"相同平行读数的反向份数"变量上做
  带度数上下界剪枝的回溯，再对完整定向枚举欧拉回路），不是只尝试原批/整批
  翻转两种；
- 结果仍按**循环移位 + 整条反向互补**归并规范类别：整批翻转只是同一闭环的
  反向互补，多种等价定向**不构成歧义**；
- 响应额外含 `orientation_mode: "unknown"`，且见证逐项返回：
  - `read_orientations`：每条读数的 `read`、`original_sequence`、
    采用的 `orientation` 与 `oriented_sequence`（定向后序列）；
  - `evidence[]` 每项追加该读数的 `orientation` 与 `oriented_sequence`，
    相邻（含首尾）`overlap` 按定向后序列核对；
- 存在两个规范类别时返回**两份稳定见证**（`witnesses`，跨次运行确定）；
- 若**任何逐条定向都无法闭环**，返回 `no_solution` 且原因为
  `orientation_constraints_unsatisfiable`（附带评估的定向组数、平衡/连通
  失败计数与可复核的必要条件），**不会**沿用固定方向下的度数失衡误导分析员；
  另有定向无关的预检 `orientation_fragmented_graph`——把每个 k-mer 与其
  反向互补折叠为同一节点后读数端点仍不连通（该断裂与任何定向无关）。

## 运行

```bash
# 启动 API（端口可用 API_PORT 覆盖）
API_PORT=9090 docker compose up --build api

# 一次性 verify：等 API 健康检查通过后，执行构建自检 + 单元测试 +
# 唯一/歧义/无解三类冒烟 + 方向未知模式（启用/不启用）业务冒烟，
# 随后自行退出，退出码 0 表示全部通过
docker compose up --build --exit-code-from verify verify
echo "verify exit code: $?"
```

`docker compose up --abort-on-container-exit verify` 或查看退出码：

```bash
docker compose run verify; echo "exit=$?"
```

本地不使用容器时：

```bash
python3 app/main.py            # 启 API（API_PORT 环境变量可配）
python3 -m unittest discover -s tests -v
python3 verify.py              # 需先启动 API，API_HOST/API_PORT 可配
```

## 目录

```
app/barcode.py    核心算法：校验、建图、欧拉回路枚举、规范归一化
app/main.py       HTTP 服务（/health、/assemble）
tests/            单元测试（含全排列暴力参照交叉验证）
verify.py         一次性 verify 服务（等待就绪→构建→测试→三类冒烟→退出码）
Dockerfile        单一镜像，默认启动 API
docker-compose.yml api + 一次性 verify
```
