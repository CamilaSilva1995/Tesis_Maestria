"""Genera Articulo_RevistaCientifica_ES.docx a partir de articulo_ES.md y ajusta el
formato a las normas de la Revista Científica (Universidad Distrital):
Times New Roman 12, interlineado 1.5, A4, márgenes 2.5 cm arriba/abajo y 3 cm laterales.
Uso: python3 build_docx.py [entrada.md] [salida.docx]
"""
import re, shutil, subprocess, zipfile, os, sys

SRC = sys.argv[1] if len(sys.argv) > 1 else "articulo_ES.md"
OUT = sys.argv[2] if len(sys.argv) > 2 else "Articulo_RevistaCientifica_ES.docx"
TMP = "_pandoc_tmp.docx"

subprocess.run(["pandoc", SRC, "-o", TMP, "--from", "markdown+raw_attribute"], check=True)

FONT = 'w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="Times New Roman" w:cs="Times New Roman"'
RPR_BASE = f'<w:rFonts {FONT}/><w:sz w:val="24"/><w:szCs w:val="24"/>'

def set_style(xml, sid, ppr, rpr):
    """Reemplaza pPr y rPr de un estilo de párrafo respetando el orden del esquema."""
    pat = re.compile(r'(<w:style [^>]*w:styleId="%s"[^>]*>)(.*?)(</w:style>)' % sid, re.S)
    m = pat.search(xml)
    if not m:
        return xml
    body = re.sub(r'<w:pPr>.*?</w:pPr>', '', m.group(2), flags=re.S)
    body = re.sub(r'<w:rPr>.*?</w:rPr>', '', body, flags=re.S)
    body = re.sub(r'<w:rPr/>|<w:pPr/>', '', body)
    return xml[:m.start()] + m.group(1) + body + ppr + rpr + m.group(3) + xml[m.end():]

def fix_styles(xml):
    xml = re.sub(r'<w:rFonts [^/]*/>', f'<w:rFonts {FONT}/>', xml)
    xml = re.sub(r'<w:sz w:val="\d+"/>', '<w:sz w:val="24"/>', xml)
    xml = re.sub(r'<w:szCs w:val="\d+"/>', '<w:szCs w:val="24"/>', xml)
    body_ppr = '<w:pPr><w:spacing w:before="0" w:after="120" w:line="360" w:lineRule="auto"/><w:jc w:val="both"/></w:pPr>'
    body_rpr = f'<w:rPr>{RPR_BASE}</w:rPr>'
    for sid in ("Normal", "BodyText", "FirstParagraph", "Compact"):
        xml = set_style(xml, sid, body_ppr, body_rpr)
    head_ppr = '<w:pPr><w:keepNext/><w:spacing w:before="240" w:after="120" w:line="360" w:lineRule="auto"/><w:jc w:val="left"/></w:pPr>'
    head_rpr = f'<w:rPr>{RPR_BASE}<w:b/><w:bCs/></w:rPr>'
    for sid in ("Heading1", "Heading2", "Heading3", "Title", "Subtitle"):
        xml = set_style(xml, sid, head_ppr, head_rpr)
    cap_ppr = '<w:pPr><w:spacing w:before="120" w:after="120" w:line="360" w:lineRule="auto"/><w:jc w:val="left"/></w:pPr>'
    for sid in ("ImageCaption", "Caption", "Figure", "CaptionedFigure", "TableCaption"):
        xml = set_style(xml, sid, cap_ppr, body_rpr)
    return xml

SECT = ('<w:sectPr><w:pgSz w:w="11906" w:h="16838" w:orient="portrait"/>'
        '<w:pgMar w:top="1417" w:right="1701" w:bottom="1417" w:left="1701" w:header="708" w:footer="708" w:gutter="0"/>'
        '<w:cols w:space="708"/></w:sectPr>')

def fix_document(xml):
    xml = re.sub(r'<w:sectPr>.*?</w:sectPr>|<w:sectPr\s*/>', '', xml, flags=re.S)
    xml = xml.replace('</w:body>', SECT + '</w:body>')
    # Párrafos que contienen una imagen: centrados (spacing antes de jc, según el esquema)
    xml = re.sub(r'<w:p>(<w:pPr>(?:(?!</w:pPr>).)*</w:pPr>)?(<w:r>(?:(?!</w:r>).)*<w:drawing>)',
                 lambda m: '<w:p><w:pPr><w:spacing w:before="120" w:after="60"/><w:jc w:val="center"/></w:pPr>' + m.group(2),
                 xml, flags=re.S)
    return xml

def fix_content_types(xml):
    for ext, ct in (("jpg", "image/jpeg"), ("jpeg", "image/jpeg"), ("png", "image/png")):
        if f'Extension="{ext}"' not in xml:
            xml = xml.replace('<Default Extension="xml"', f'<Default Extension="{ext}" ContentType="{ct}"/><Default Extension="xml"', 1)
    return xml

zin = zipfile.ZipFile(TMP)
zout = zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED)
for item in zin.infolist():
    data = zin.read(item.filename)
    if item.filename == "word/styles.xml":
        data = fix_styles(data.decode("utf8")).encode("utf8")
    elif item.filename == "word/document.xml":
        data = fix_document(data.decode("utf8")).encode("utf8")
    elif item.filename == "[Content_Types].xml":
        data = fix_content_types(data.decode("utf8")).encode("utf8")
    zout.writestr(item, data)
zout.close(); zin.close(); os.remove(TMP)
print("OK", OUT)
