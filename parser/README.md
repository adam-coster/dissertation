# Dissertation LaTeX Parser

This folder contains Node.js code for parsing the various LaTeX and attendant files making up the source of my dissertation and converting them into a more portable data format, for use in creating alternate versions of the dissertation (e.g. for HTML).

This code is VERY specific to the details of this particular LaTeX project, and gets to make assumptions that wouldn't be *generally* true in LaTeX projects.

## Usage

Run scripts from the repo root (e.g. `node parser/parse.glossary.mjs`). Requires Node.js 20+.

## Glossary Parser

