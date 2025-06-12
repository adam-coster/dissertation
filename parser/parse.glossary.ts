import { ok } from "node:assert";
import fsp from "node:fs/promises";
import { parseInline } from "./parse.utility.js";
import type { GlossaryEntry } from "./types.js";

const glossaryText = await fsp.readFile("../glossary.tex", "utf-8");

// Global split, getting rid of preceding whitespace including newlines
const rawEntries = glossaryText
  .split(/[\s\r\n]*\\newglossaryentry/g)
  .filter(Boolean);

console.log(rawEntries);

for (const rawEntry of rawEntries) {
  let { id, rawBody } = rawEntry.match(/\{(?<id>[^}]+)\}(?<rawBody>.*)/s)!
    .groups as { id: string; rawBody: string };
  ok(id, `Glossary entry ID is required: ${rawEntry}`);
  ok(rawBody, `Glossary entry body is required: ${rawEntry}`);
  let [rawName, ...rawDescriptionLines] = rawBody
    .trim()
    .replace(/^\{\s*(.*)\s*\}$/s, "$1")
    .split(/[\r\n]+/);
  rawName = rawName
    .trim()
    .match(/^name\s*=\s*(.*?)\s*$/)![1]
    .trim();
  const { itStart, name } = rawName.match(
    /^(?<itStart>\\textit\{)?(?<name>.*?)(?<itEnd>\})?\s*,\s*$/
  )!.groups as {
    itStart?: string;
    name: string;
    itEnd?: string;
  };
  const entry: GlossaryEntry = {
    type: "glossary-entry",
    id,
    name,
    italic: !!itStart,
    description: [],
  };

  // Pull out the plural form if it exists. All are "plural=\textit{...}"
  const pluralLine = rawDescriptionLines.findIndex((line) =>
    line.trim().startsWith("plural=")
  );
  if (pluralLine !== -1) {
    const pluralForm = rawDescriptionLines[pluralLine]
      .replace(/^\s*plural\s*=\s*\\textit\{(.*?)\}\s*$/, "$1")
      .trim();
    entry.plural = {
      type: "italic",
      children: [{ type: "text", value: pluralForm }],
    };
    rawDescriptionLines.splice(pluralLine, 1);
  }

  const rawDescription = rawDescriptionLines
    .map((l) => l.trim())
    .join(" ")
    .replace(/^description\s*=\s*\{(.*?)\s*\}(\s*,)?\s*$/, "$1")
    .trim();
  // Split up the description into nodes, since it can contain italics and references
  entry.description = parseInline(rawDescription);

  if (entry.plural) {
    console.log(entry);
  }
}
