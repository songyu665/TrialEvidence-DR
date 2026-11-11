import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def load(name: str) -> dict:
    return json.loads(project_path("reports", name).read_text(encoding="utf-8"))


def main() -> None:
    trials = load("e0_trial_sample_audit.json")
    negative = load("e0_aggregated_n_strict_candidates.json")
    review = load("e0_conservative_rule_review.json")
    positive = load("e0_chembl_approved_candidates.json")
    conflicts = load("e0_conflict_and_coverage.json")
    gate = load("e0_rule_based_gate.json")
    completed = load("e0_completed_results_audit.json")
    review_pool = load("e0_completed_n_review_mapping.json")
    ai_review = load("e0_ai_assisted_review_200_audited.json")
    integrity = load("raw_download_integrity.json")
    counts = completed["audit_counts"]
    lines = [
        "# TrialEvidence-DR E0 数据可行性审计",
        "",
        f"生成时间：{datetime.now(timezone.utc).isoformat()}",
        "",
        "> 本报告区分文本候选、规则构建标签和领域确认标签。当前没有临床金标准标签，也没有模型性能结论。",
        "",
        "## 数据完整性与范围",
        "",
        "| 数据集合 | 当前记录数 | 完整下载 | 校验 |",
        "|---|---:|---|---|",
    ]
    for name, result in integrity["datasets"].items():
        lines.append(f"| {name} | {result['records_verified']} | {result['download_complete']} | {result['status']} |")
    lines += [
        "",
        "## 终止试验阴性证据",
        "",
        f"- 明确查询集合：{trials['records_audited']}条，范围为{trials['sample_design']}。",
        f"- 明确 futility/缺乏疗效文本候选：{trials['candidate_buckets'].get('N_strict_candidate_requires_three_layer_review', 0)}个试验。",
        f"- 可追溯映射并聚合后的唯一候选配对：{negative['unique_pair_candidates']}，覆盖{negative['unique_drugs']}个药物和{negative['unique_diseases']}个疾病节点。",
        f"- 保守干预归属和上下文规则通过：{review['passed_trial_pair_records']}条试验—配对记录，{review['passed_unique_pair_candidates']}个唯一配对。",
        "- 通过者是规则构建候选，不是领域专家确认的临床阴性。",
        "",
        "## 批准阳性候选",
        "",
        f"- ChEMBL 原始记录：{positive['source_records']}。",
        f"- 最高阶段4且具有 DailyMed/EMA/FDA 引用的去重候选：{positive['approved_candidate_pairs']}。",
        f"- 与阴性候选重叠：{conflicts['positive_negative_conflicting_pairs']}对；重叠需要上下文复核，不能自动判定某一方错误。",
        "",
        "## 已完成且有结果的试验",
        "",
        f"- 当前审计记录：{counts.get('studies', 0)}/{completed['source_total_reported_by_api']}，完整下载={completed['source_download_is_complete']}。",
        f"- 有主要终点分析：{counts.get('with_primary_analysis', 0)}；有 P 值：{counts.get('with_primary_pvalue', 0)}；有置信区间：{counts.get('with_primary_confidence_interval', 0)}。",
        f"- 含非劣效分析：{counts.get('with_noninferiority_analysis', 0)}。",
        "- 所有记录仅进入 N_review；P>0.05 不等于无效，未自动增加 N_strict。",
        f"- 可映射的唯一 N_review 配对：{review_pool['mapped_unique_n_review_pairs']}，覆盖{review_pool['mapped_unique_drugs']}个药物和{review_pool['mapped_unique_diseases']}个疾病。",
        "",
        "## AI 辅助分层复核",
        "",
        f"- 高文本优先级：0/{ai_review['strata']['high']['n']}支持阴性候选；该层为全量纳入。",
        f"- 随机中等优先级：{ai_review['strata']['medium']['decisions']['supports_negative_candidate']}/{ai_review['strata']['medium']['n']}支持，比例{ai_review['strata']['medium']['support_fraction']:.2%}，Wilson 95%区间 {ai_review['strata']['medium']['wilson_95'][0]:.2%}—{ai_review['strata']['medium']['wilson_95'][1]:.2%}。",
        "- 两层不能合并为总体估计；父级审查已将一个联合/加用方案纠正为不适用于单药标签；AI 筛查不等于医学领域确认。",
        "",
        "## 当前门控结论",
        "",
        "**E0 未通过，不进入完整模型比较。**",
        "",
        f"1. 保守规则后有{review['passed_unique_pair_candidates']}个唯一候选；其中{gate['positive_overlapping_pairs']}个与批准阳性重叠，剩余{gate['nonoverlapping_rule_based_pairs']}个非重叠候选，低于项目管理门槛200。", 
        "2. 三层文本与规则审查不能替代医学/药理领域复核。",
        "3. 历史图谱、历史特征和信息首次公开时间尚未完成可恢复性审计。",
        "4. 已完成结果试验可能扩展 N_review，但需要效应方向、对照、估计目标和多重性裁决。",
        "",
        "## 允许的下一步",
        "",
        "- 完成已完成结果试验的全量审计和保守裁决，不降低阴性语义标准凑数量。",
        "- 若可靠阴性仍不足，转向小规模标签资源、标签质量和选择偏倚研究，或重新定义候选总体。",
        "- 只有在标签覆盖、复核、冲突和划分协议均冻结后，才开展 PN/nnPU/PNU/PUbN 最小对照。",
        "",
    ]
    output = project_path("reports", "E0_数据可行性审计.md")
    output.write_text("\n".join(lines), encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
