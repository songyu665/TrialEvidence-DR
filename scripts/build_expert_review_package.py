import csv
import json
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path

OUT = project_path("expert_review")
REVIEW_FIELDS = ["reviewer", "review_date", "extraction_fidelity", "intervention_attribution", "context_aggregation", "statistical_semantics", "clinical_semantics", "final_adjudication", "exclusion_reason", "review_notes"]


def read_csv(name: str) -> list[dict]:
    with project_path("data", "interim", name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name: str, rows: list[dict]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    counts = {}
    rule_rows = read_csv("n_strict_rule_based_pairs.csv")
    output = []
    for index, row in enumerate(rule_rows, 1):
        output.append({
            "review_id": f"N_RULE_{index:04d}", "drug_node_name": row["drug_node_name"],
            "disease_node_name": row["disease_node_name"], "trial_ids": row["trial_ids"],
            "raw_evidence": row["raw_evidence"], "evidence_urls": row["evidence_urls"],
            "mapping_transforms": row["drug_mapping_transforms"], "positive_overlap": row["positive_overlap"],
            **{field: "" for field in REVIEW_FIELDS},
        })
    write_csv("01_rule_based_negative_pairs.csv", output)
    counts["rule_based_negative_pairs"] = len(output)

    conflicts = read_csv("positive_negative_conflict_review.csv")
    output = []
    for index, row in enumerate(conflicts, 1):
        output.append({
            "review_id": f"CONFLICT_{index:04d}", "drug_node_name": row["drug_node_name"],
            "disease_node_name": row["disease_node_name"], "negative_trial_ids": row["negative_trial_ids"],
            "negative_raw_evidence": row["negative_raw_evidence"],
            "positive_reference_urls": row["positive_reference_urls"],
            **{field: "" for field in REVIEW_FIELDS},
        })
    write_csv("02_positive_negative_conflicts.csv", output)
    counts["positive_negative_conflicts"] = len(output)

    audited = read_csv("e0_ai_assisted_review_200_audited.csv")
    supported = [row for row in audited if row["ai_screening_decision"] == "supports_negative_candidate"]
    output = []
    for index, row in enumerate(supported, 1):
        output.append({
            "review_id": f"AI_SUPPORT_{index:04d}", "pair_id": row["pair_id"],
            "screening_reason": row["reason"], "evidence_nct_ids": row["evidence_nct_ids"],
            "parent_audit_note": row["parent_audit_note"],
            **{field: "" for field in REVIEW_FIELDS},
        })
    write_csv("03_ai_supported_cases_for_expert_review.csv", output)
    counts["ai_supported_cases"] = len(output)

    dates = {row["pair_id"]: row for row in read_csv("dailymed_label_date_audit.csv")}
    events = read_csv("post_cutoff_positive_event_candidates.csv")
    output = []
    for index, row in enumerate(events, 1):
        pair_id = f"{row['drug_node_index']}|{row['disease_node_index']}"
        date_info = dates.get(pair_id, {})
        output.append({
            "review_id": f"EVENT_{index:04d}", "pair_id": pair_id,
            "drug_node_name": row["drug_node_name"], "disease_node_name": row["disease_node_name"],
            "event_category": row["event_category"], "new_drugind_ids": row["new_drugind_ids_relative_to_chembl33"],
            "existing_drugind_ids": row["existing_drugind_ids_in_chembl33"],
            "first_approval_years": row["first_approval_years"], "reference_urls": row["reference_urls"],
            "dailymed_label_dates": date_info.get("label_published_dates", ""),
            "dailymed_date_status": date_info.get("date_audit_status", "not_available"),
            "original_indication_event_date": "", "event_date_source": "", "event_type": "",
            "reviewer": "", "review_date": "", "review_notes": "",
        })
    write_csv("04_post_cutoff_event_date_review.csv", output)
    counts["post_cutoff_event_candidates"] = len(output)

    instructions = """# 外部专家复核说明

本包不包含临床金标准。请独立判断每项记录，不要因为已有 AI 或规则结论而默认接受。

## 复核层级

1. 抽取忠实性：原始记录是否确实表达 futility、缺乏疗效或相应事件。
2. 干预归属：结果是否归属于目标药物、目标治疗臂和当前方案，而非其他臂或共同背景治疗。
3. 上下文聚合：疾病亚型、人群、治疗线次、单药/联合方案是否支持当前配对。
4. 统计语义：效应方向、对照、估计目标、区间、效能和多重性是否支持候选解释。P>0.05 单独不等于无效。
5. 临床语义：能否作为当前研究定义下的可靠阴性或阳性；不能确认时保留为 U。

## 文件

- 01：36个规则构建阴性候选，其中2个与批准阳性重叠；均需独立医学/药理复核。
- 02：批准阳性与试验阴性重叠，需判断是否为不同适应证上下文或真正冲突。
- 03：AI 辅助筛查支持的个案，不能视为已确认。
- 04：潜在后续事件。ChEMBL 新记录 ID 和药物首次批准年份不等于适应证批准日期。DailyMed published_date 只是标签发布日期。

请填写复核者、日期、各层判断、最终裁决、排除理由和备注。不要根据预期研究结论反向选择标签。
"""
    (OUT / "README.md").write_text(instructions, encoding="utf-8")
    manifest = {"created_at": datetime.now(timezone.utc).isoformat(), "counts": counts, "status": "pending_independent_domain_review"}
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    zip_path = project_path("TrialEvidence_DR_expert_review.zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.iterdir()):
            archive.write(path, arcname=f"expert_review/{path.name}")
    with zipfile.ZipFile(zip_path) as archive:
        if archive.testzip():
            raise ValueError("Expert-review ZIP failed integrity check")
    print(json.dumps({"zip": str(zip_path), "counts": counts, "validation": "passed"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
