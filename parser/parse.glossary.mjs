import fsp from "node:fs/promises";

const glossaryText = await fsp.readFile("./glossary.tex", "utf-8");
