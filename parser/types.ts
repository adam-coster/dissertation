export interface ParagraphNode {
  type: "paragraph";
  value: (TextNode | ItalicNode | GlossaryReference)[];
}

export interface TextNode {
  type: "text";
  value: string;
}

export interface ItalicNode {
  type: "italic";
  value: TextNode;
}

export interface GlossaryEntry {
  type: "glossary-entry";
  id: string;
  name: string;
  italic?: boolean;
  description: (TextNode | ItalicNode | GlossaryReference)[];
}

export interface GlossaryReference {
  type: "glossary-reference";
  id: string;
}
