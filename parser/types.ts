export interface ParagraphNode {
  type: "paragraph";
  value: Inline[];
}

export interface Text {
  type: "text";
  value: string;
}

export interface LeftQuote {
  type: "left-quote";
  double?: boolean;
}

export interface RightQuote {
  type: "right-quote";
  double?: boolean;
}

export type Inline =
  | Text
  | ItalicNode
  | GlossaryReference
  | BoldNode
  | LeftQuote
  | RightQuote;

export interface ItalicNode {
  type: "italic";
  children: Inline[];
}

export interface BoldNode {
  type: "bold";
  children: Inline[];
}

export interface GlossaryEntry {
  type: "glossary-entry";
  id: string;
  name: string;
  plural?: ItalicNode;
  italic?: boolean;
  description: Inline[];
}

export interface GlossaryReference {
  type: "glossary-reference";
  id: string;
}
