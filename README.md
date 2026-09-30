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

**方向未知批次**（`POST /assemble` 加 `"unknown_orientation": true`）：
部分测序批次会丢失单条读数的链方向。启用后每条读数可独立按**原序列或其
反向互补**参与拼接，但每个输入序号仍只能使用一次；"逐条定向"与"闭环
次序"在同一次 DFS 中**联合求解**（不是只尝试整批翻转）。无解时返回方向
松弛后的可复核原因，而非沿用固定方向下会误导分析员的度数失衡：

- `orientation_group_imbalance`：把 k-mer 与它的反向互补归为同一配对组后，
  存在关联读数条数为奇数的组（翻转只改变组间方向、不改变组，奇度无法消除）；
- `orientation_group_fragmented`：配对组落在多个互不连通的分量中；
- `orientation_constraints_unsatisfiable`：通过上述必要条件但联合穷尽枚举
  仍无法闭环，附搜索统计与 `dead_end_examples`（停在哪个 k-mer、剩余读数
  各自两种定向的候选边），可逐项复核。

规范等价归一化：

1. 环的循环移位；
2. 整条环的反向互补；
3. 相同读数（平行边）互换不产生新类别（枚举时同标签只取当前最小输入序号，
   因而**使用次序按输入序号稳定分配**且跨次运行确定）。

枚举从固定起点出发做带桥剪枝的欧拉回路 DFS，按规范代表去重，找到 2 个
不同等价类即停。约束（6~24 条、长度 3~8）下毫秒级完成。

## API

零第三方依赖（Python 标准库）。

- `GET /health` → `{"status":"ok",...}`
- `POST /assemble`，请求体：

```json
{ "sequences": ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"] }
```

方向未知批次追加布尔字段（缺省/false 时与历史版本完全兼容）：

```json
{
  "sequences": ["ATT", "ATC", "TCG", "CGC", "GCA", "TTG"],
  "unknown_orientation": true
}
```

响应 `status` 为 `unique` / `ambiguous` / `no_solution`（非法输入返回 400，
非法 `unknown_orientation` 值同样返回 400）。

`unique` 响应包含：

- `canonical_barcode`：规范条码（环旋转+反向互补下的字典序最小代表）；
- `barcode`：与 `order` 同方向的环；
- `order`：使用次序，元素为 **1 起始的输入序号**（每个序号恰好出现一次）；
- `evidence`：相邻重叠证据（`prev`/`next` 序号、`overlap`、`overlap_length = L-1`、
  `appended_base`），最后一条与首条首尾相接。

`ambiguous` 响应在 `witnesses` 中给出两条不同规范见证（结构同上）。

方向未知模式额外包含 `orientation_mode: "unknown"` 以及：

- `orientations`：与 `order` 对齐的采用方向（`forward` / `reverse_complement`）；
- 每条证据额外给出 `orientation`、`original_sequence`、`oriented_sequence`，
  即逐项返回**采用方向、定向后序列及相邻重叠**；`barcode` 已对齐到规范代表。

规范归并仍按循环移位与整条反向互补进行：同一条码的多种等价定向（如周期
环、回文读数）已合并，不会制造歧义；存在两个条码类别时 `witnesses` 恒为
两份稳定见证（重复调用次序确定）。

## 运行

```bash
# 启动 API（端口可用 API_PORT 覆盖）
API_PORT=9090 docker compose up --build api

# 一次性 verify：等 API 健康检查通过后，执行构建自检 + 单元测试 +
# 唯一/歧义/无解三类冒烟，随后自行退出，退出码 0 表示全部通过
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
