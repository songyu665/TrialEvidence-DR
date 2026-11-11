import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from trialevidence.paths import project_path


def clean(text: str) -> str:
    text = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r"\1 (\2)", text)
    return text.replace("**", "").replace("`", "").strip()


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instruction = OxmlElement("w:instrText")
    instruction.set(qn("xml:space"), "preserve")
    instruction.text = "PAGE"
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend([begin, instruction, end])


def main() -> None:
    source = project_path("manuscript", "TrialEvidence_DR_manuscript.md")
    output = project_path("manuscript", "TrialEvidence_DR_manuscript.docx")
    document = Document()
    section = document.sections[0]
    section.top_margin = section.bottom_margin = Inches(0.85)
    section.left_margin = section.right_margin = Inches(0.9)
    normal = document.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(11)
    normal.paragraph_format.line_spacing = 1.15
    normal.paragraph_format.space_after = Pt(6)
    add_page_number(section.footer.paragraphs[0])

    for raw in source.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        heading = re.match(r"^(#{1,6})\s+(.*)", line)
        if heading:
            level = len(heading.group(1))
            text = clean(heading.group(2))
            document.add_heading(text, level=0 if level == 1 else min(level - 1, 4))
            continue
        if line.startswith("- "):
            document.add_paragraph(clean(line[2:]), style="List Bullet")
            continue
        image_paths = re.findall(r"(?:reports|data)/[^\s`,]+?\.png", line)
        if image_paths:
            for relative in image_paths:
                path = project_path(relative)
                document.add_picture(str(path), width=Inches(6.4))
                document.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
            caption = clean(line)
            if caption:
                document.add_paragraph(caption, style="Caption")
            continue
        document.add_paragraph(clean(line))

    temporary = output.with_suffix(".docx.tmp")
    document.save(temporary)
    temporary.replace(output)
    verified = Document(output)
    report = {
        "output": str(output), "paragraphs": len(verified.paragraphs),
        "embedded_figures": len(verified.inline_shapes), "bytes": output.stat().st_size,
        "title": verified.paragraphs[0].text,
        "last_paragraph": verified.paragraphs[-1].text,
    }
    if not report["embedded_figures"] or report["paragraphs"] < 50:
        raise ValueError("DOCX verification failed")
    print(report)


if __name__ == "__main__":
    main()
