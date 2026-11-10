# Preserve query rows when adding sample metadata

Inverted queries retain one row per query: `any-bins` and `all-bins` list comma-separated reference matches, and `match-count` retains its numeric matrix columns. Append optional metadata to names with a semicolon, including reference names in matrix headers, and reject tabs, commas, and semicolons in metadata so consumers can split fields without escaping. This preserves the existing row structure; `--suppress-metadata` provides unannotated names for consumers that require them.
