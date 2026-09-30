"""一次性 verify 容器：

  1. 等待 API /health 就绪（内部轮询，满足"verify 服务等待 API 就绪"）；
  2. 构建自检（compileall 字节码编译）；
  3. 运行单元测试；
  4. 通过 HTTP 对 /assemble 做唯一、歧义、无解三类冒烟；
  5. 对方向未知（unknown_orientation）模式做新旧兼容、唯一恢复、整批翻转
     同类、歧义双见证、不可满足原因与折叠图断裂等业务冒烟；
  6. 汇总结果后自行退出，全部通过退出码 0，否则非零。
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request

API_HOST = os.environ.get("API_HOST", "api")
API_PORT = os.environ.get("API_PORT", "8080")
BASE = f"http://{API_HOST}:{API_PORT}"
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "app"))

failures: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    mark = "PASS" if ok else "FAIL"
    print(f"[{mark}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        failures.append(name)
    return ok


def wait_ready(timeout: float = 60.0) -> bool:
    deadline = time.time() + timeout
    last = ""
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(f"{BASE}/health", timeout=3) as resp:
                if resp.status == 200:
                    body = json.loads(resp.read().decode())
                    return body.get("status") == "ok"
        except (urllib.error.URLError, OSError, json.JSONDecodeError) as exc:
            last = str(exc)
            time.sleep(0.5)
    print(f"等待 API 就绪超时: {last}")
    return False


def post_assemble(
    reads: list[str], unknown_orientation: bool | None = None
) -> tuple[int, dict]:
    payload: dict = {"sequences": reads}
    if unknown_orientation is not None:
        payload["unknown_orientation"] = unknown_orientation
    req = urllib.request.Request(
        f"{BASE}/assemble",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            return resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode())


def smoke_unique() -> None:
    reads = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]
    code, body = post_assemble(reads)
    ok = (
        code == 200
        and body.get("status") == "unique"
        and body.get("canonical_barcode") == "AATCGC"
        and sorted(body.get("order", [])) == list(range(1, 7))
        and len(body.get("evidence", [])) == 6
        and all(e["overlap_length"] == 2 for e in body["evidence"])
    )
    check("冒烟-唯一拼接", bool(ok),
          f"status={body.get('status')} barcode={body.get('canonical_barcode')}")

    # 反向互补整组提交必须得到同一规范条码（等价归一化）
    from barcode import reverse_complement
    code2, body2 = post_assemble([reverse_complement(s) for s in reads])
    check("冒烟-反向互补等价",
          code2 == 200 and body2.get("canonical_barcode") == "AATCGC",
          f"barcode={body2.get('canonical_barcode')}")


def smoke_ambiguous() -> None:
    reads = ["AAC", "ACC", "AAG", "AGC", "GCC",
             "CCA", "CAA", "CCG", "CGA", "GAA"]
    code, body = post_assemble(reads)
    witnesses = body.get("witnesses", [])
    barcodes = {w.get("canonical_barcode") for w in witnesses}
    ok = (
        code == 200
        and body.get("status") == "ambiguous"
        and len(witnesses) == 2
        and len(barcodes) == 2
        and all(sorted(w.get("order", [])) == list(range(1, 11)) for w in witnesses)
        and all(len(w.get("evidence", [])) == 10 for w in witnesses)
    )
    check("冒烟-歧义双见证", bool(ok),
          f"witnesses={sorted(barcodes)}")


def smoke_no_solution() -> None:
    # 两个原因同时存在：度数失衡 + 多分量
    reads = ["AAA", "AAA", "AAC", "CCC", "GGG", "TTT"]
    code, body = post_assemble(reads)
    reasons = {r.get("code") for r in body.get("reasons", [])}
    ok = (
        code == 200
        and body.get("status") == "no_solution"
        and {"degree_imbalance", "fragmented_graph"} <= reasons
    )
    check("冒烟-无解(度数失衡+非零片段不连通)", bool(ok),
          f"reasons={sorted(reasons)}")

    # 仅多分量但各自平衡
    reads2 = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
    code2, body2 = post_assemble(reads2)
    reasons2 = {r.get("code") for r in body2.get("reasons", [])}
    check("冒烟-无解(多分量各自平衡)",
          code2 == 200 and body2.get("status") == "no_solution"
          and reasons2 == {"fragmented_graph"},
          f"reasons={sorted(reasons2)}")


def smoke_invalid() -> None:
    code, body = post_assemble(["AAA"] * 5)  # 少于 6 条
    check("冒烟-非法输入 400", code == 400 and body.get("status") == "invalid_input",
          f"code={code}")


def smoke_unknown_orientation() -> None:
    from barcode import reverse_complement

    reads = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]

    # 1) 未启用标志：响应结构与旧模式完全兼容（无方向字段）
    code, body = post_assemble(reads)
    check("冒烟-旧模式默认结构兼容",
          code == 200 and "orientation_mode" not in body
          and "orientation" not in body["evidence"][0],
          f"status={body.get('status')}")

    # 2) 部分读数丢失链方向（随机 RC），联合定向后恢复唯一规范条码
    flipped = [reverse_complement(s) for s in
               [reads[0], reads[2], reads[4]]]
    batch = [flipped[0], reads[1], flipped[1], reads[3], flipped[2], reads[5]]
    code, body = post_assemble(batch, unknown_orientation=True)
    ok = (
        code == 200
        and body.get("status") == "unique"
        and body.get("orientation_mode") == "unknown"
        and body.get("canonical_barcode") == "AATCGC"
        and sorted(body.get("order", [])) == list(range(1, 7))
        and len(body.get("read_orientations", [])) == 6
        and len(body.get("evidence", [])) == 6
    )
    # 逐项见证：方向、定向后序列、相邻重叠闭环一致
    oriented = {
        ro["read"]: ro["oriented_sequence"]
        for ro in body.get("read_orientations", [])
    }
    for ev in body.get("evidence", []):
        a, b = oriented.get(ev["prev"]), oriented.get(ev["next"])
        ok = ok and a is not None and b is not None
        ok = ok and ev.get("oriented_sequence") == a
        ok = ok and ev["overlap"] == a[1:] == b[:-1]
        ok = ok and ev.get("orientation") in ("forward", "reverse_complement")
    expected_dirs = {1: "reverse_complement", 2: "forward",
                     3: "reverse_complement", 4: "forward",
                     5: "reverse_complement", 6: "forward"}
    got_dirs = {ro["read"]: ro["orientation"]
                for ro in body.get("read_orientations", [])}
    check("冒烟-方向未知唯一恢复(逐项定向见证)",
          bool(ok) and got_dirs == expected_dirs,
          f"barcode={body.get('canonical_barcode')} dirs={got_dirs}")

    # 3) 整批反向互补：同一规范类，不构成歧义
    code, body = post_assemble(
        [reverse_complement(s) for s in reads], unknown_orientation=True
    )
    check("冒烟-方向未知整批翻转同类",
          code == 200 and body.get("status") == "unique"
          and body.get("canonical_barcode") == "AATCGC",
          f"status={body.get('status')}")

    # 4) 两个规范类别：返回两份稳定见证
    theta = ["AAC", "ACC", "AAG", "AGC", "GCC",
             "CCA", "CAA", "CCG", "CGA", "GAA"]
    code, body = post_assemble(theta, unknown_orientation=True)
    witnesses = body.get("witnesses", [])
    barcodes = sorted(w.get("canonical_barcode") for w in witnesses)
    ok = (
        code == 200 and body.get("status") == "ambiguous"
        and len(witnesses) == 2 and len(set(barcodes)) == 2
        and all(sorted(w.get("order", [])) == list(range(1, 11))
                for w in witnesses)
        and all(len(w.get("read_orientations", [])) == 10 for w in witnesses)
    )
    code2, body2 = post_assemble(theta, unknown_orientation=True)
    stable = [w.get("canonical_barcode") for w in body2.get("witnesses", [])]
    check("冒烟-方向未知歧义双见证且稳定",
          bool(ok) and sorted(stable) == barcodes,
          f"witnesses={barcodes}")

    # 5) 任何逐条定向都无法闭环：方向约束不可满足，而非度数失衡
    unsat = ["AAA", "AAA", "AAC", "ACA", "CAA", "CAA"]
    code, body = post_assemble(unsat, unknown_orientation=True)
    codes = {r.get("code") for r in body.get("reasons", [])}
    check("冒烟-方向未知不可满足(非度数失衡)",
          code == 200 and body.get("status") == "no_solution"
          and codes == {"orientation_constraints_unsatisfiable"},
          f"reasons={sorted(codes)}")

    # 6) 折叠图不连通（定向无关的预检）
    frag = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
    code, body = post_assemble(frag, unknown_orientation=True)
    codes = {r.get("code") for r in body.get("reasons", [])}
    check("冒烟-方向未知折叠图不连通",
          code == 200 and body.get("status") == "no_solution"
          and codes == {"orientation_fragmented_graph"},
          f"reasons={sorted(codes)}")

    # 7) 非法标志类型返回 400
    req = urllib.request.Request(
        f"{BASE}/assemble",
        data=json.dumps({"sequences": reads,
                         "unknown_orientation": "yes"}).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            bad_code = resp.status
    except urllib.error.HTTPError as exc:
        bad_code = exc.code
    check("冒烟-方向标志非法 400", bad_code == 400, f"code={bad_code}")


def main() -> int:
    print(f"verify: target={BASE}", flush=True)
    if not check("等待 API 就绪", wait_ready()):
        return 1

    build = subprocess.run(
        [sys.executable, "-m", "compileall", "-q", "app", "tests", "verify.py"],
        cwd=ROOT,
    )
    check("构建自检 compileall", build.returncode == 0)

    tests = subprocess.run(
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
        cwd=ROOT,
    )
    check("单元测试套件", tests.returncode == 0)

    smoke_unique()
    smoke_ambiguous()
    smoke_no_solution()
    smoke_invalid()
    smoke_unknown_orientation()

    if failures:
        print(f"\nverify 失败 {len(failures)} 项: {failures}", flush=True)
        return 1
    print("\nverify 全部通过", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
