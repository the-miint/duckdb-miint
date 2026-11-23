#### `parse_gff_attributes(kvp_string)`

<!-- manual:url /reading/#parse_gff_attributeskvp_string -->
<!-- manual:description -->
Parse a GFF3 attribute string into a map of keys to values.

`read_gff` applies this macro to column 9. You can also call it directly.

**Parameters:**

- `kvp_string` (VARCHAR): A semicolon-separated `key=value` string (a GFF3 column-9 value).

**Returns:** `MAP(VARCHAR, VARCHAR)` — the parsed key/value pairs.

<!-- manual:example -->
```sql
SELECT parse_gff_attributes('ID=gene1;Name=TEST1;biotype=protein_coding');
-- {'ID': 'gene1', 'Name': 'TEST1', 'biotype': 'protein_coding'}
```

- Access values using bracket notation: `attributes['ID']`
