# Dissertation LaTeX Parser

This folder contains Node.js code for parsing the various LaTeX and attendant files making up the source of my dissertation and converting them into a more portable data format, for use in creating alternate versions of the dissertation (e.g. for HTML).

This code is VERY specific to the details of this particular LaTeX project, and gets to make assumptions that wouldn't be *generally* true in LaTeX projects.

## Usage

Run `package.json` scripts from this folder. Requires Node.js 22+ for Typescript support.

## Glossary Parser

`npm run parse:glossary`


