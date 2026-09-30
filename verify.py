"""一次性 verify 容器：

  1. 等待 API /health 就绪（内部轮询，满足"verify 服务等待 API 就绪"）；
  2. 构建自检（compileall 字节码编译）；
  3. 运行单元测试；
  4. 通过 HTTP 对 /assemble 做唯一、歧义、无解三类冒烟（固定方向模式）；
  5. 方向未知模式冒烟：逐条翻转恢复、同条码等价定向、组级不可满足、
     非法开关值 400、默认行为兼容；
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

    # 1. 逐条翻转恢复（联合定向，不是整批翻转）：4 条被翻转
    reads = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]
    flipped = [reverse_complement(s) if i in (0, 3, 4, 5) else s
               for i, s in enumerate(reads)]
    code, body = post_assemble(flipped, unknown_orientation=True)
    ev = body.get("evidence", [])
    ok = (
        code == 200
        and body.get("status") == "unique"
        and body.get("canonical_barcode") == "AATCGC"
        and sorted(body.get("order", [])) == list(range(1, 7))
        and len(body.get("orientations", [])) == 6
        and all(e.get("orientation") in ("forward", "reverse_complement")
                for e in ev)
        and all(e.get("oriented_sequence") and e.get("overlap")
                for e in ev)
    )
    check("冒烟-方向未知(逐条翻转可闭环+逐项见证)", bool(ok),
          f"status={body.get('status')} barcode={body.get('canonical_barcode')}")

    # 2. 同一条码的多种等价定向不得制造歧义（ATATAT 周期环）
    atat = ["ATA", "TAT", "ATA", "TAT", "ATA", "TAT"]
    code2, body2 = post_assemble(atat, unknown_orientation=True)
    check("冒烟-方向未知(等价定向合并,不歧义)",
          code2 == 200 and body2.get("status") == "unique"
          and body2.get("canonical_barcode") == "ATATAT",
          f"status={body2.get('status')}")

    # 3. 方向约束不可满足：给配对组奇度可复核原因，且不报固定方向度数失衡
    odd = ["AAT", "TAA", "AAT", "TAA", "AAT", "TAA"]
    code3, body3 = post_assemble(odd, unknown_orientation=True)
    reasons3 = {r.get("code") for r in body3.get("reasons", [])}
    check("冒烟-方向未知(奇度配对组不可满足)",
          code3 == 200 and body3.get("status") == "no_solution"
          and reasons3 == {"orientation_group_imbalance"}
          and "degree_imbalance" not in reasons3,
          f"reasons={sorted(reasons3)}")

    # 4. 通过预检但穷尽枚举无解：dead_end 可逐项复核
    unsat = ['GTC', 'TCT', 'ACG', 'GCG', 'GCC', 'AGG']
    code4, body4 = post_assemble(unsat, unknown_orientation=True)
    reasons4 = {r.get("code") for r in body4.get("reasons", [])}
    dead = body4.get("reasons", [{}])[0].get("dead_end_examples", [])
    check("冒烟-方向未知(穷尽无解可复核,非度数失衡)",
          code4 == 200 and body4.get("status") == "no_solution"
          and reasons4 == {"orientation_constraints_unsatisfiable"}
          and len(dead) >= 1 and all("candidate_edges" in d for d in dead),
          f"reasons={sorted(reasons4)}")

    # 5. 两类别返回两份稳定见证（即便允许逐条定向仍歧义）
    ambig = ['ATT', 'AGA', 'TTA', 'TAA', 'GAT', 'TAT', 'TAT', 'AAG']
    code5, body5 = post_assemble(ambig, unknown_orientation=True)
    witnesses = body5.get("witnesses", [])
    barcodes = {w.get("canonical_barcode") for w in witnesses}
    code5b, body5b = post_assemble(ambig, unknown_orientation=True)
    stable = [w.get("order") for w in body5b.get("witnesses", [])]
    check("冒烟-方向未知(两类别双见证且稳定)",
          code5 == 200 and body5.get("status") == "ambiguous"
          and len(witnesses) == 2 and len(barcodes) == 2
          and all(sorted(w.get("order", [])) == list(range(1, 9))
                  for w in witnesses)
          and stable == [w.get("order") for w in witnesses],
          f"witnesses={sorted(barcodes)}")

    # 6. 默认（不传）与显式 false 行为完全兼容，且 true/false 结果可不同
    c6a, b6a = post_assemble(flipped)
    c6b, b6b = post_assemble(flipped, unknown_orientation=False)
    check("冒烟-默认行为兼容(不传=false)",
          c6a == 200 and c6b == 200 and b6a == b6b
          and "orientation_mode" not in b6a,
          f"status={b6a.get('status')}")

    # 7. 非法开关值返回 400
    req = urllib.request.Request(
        f"{BASE}/assemble",
        data=json.dumps({"sequences": reads,
                         "unknown_orientation": "maybe"}).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            c7, b7 = resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        c7, b7 = exc.code, json.loads(exc.read().decode())
    check("冒烟-非法 unknown_orientation 值 400",
          c7 == 400 and b7.get("status") == "invalid_input",
          f"code={c7}")


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
