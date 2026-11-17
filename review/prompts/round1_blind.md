# Round 1: Blind Characterisation Prompt

You are given a list of sky positions with Gaia DR3 source IDs.
For each one, use public archive data to describe the object at that
position across all wavelengths you can access. Report measurements
with uncertainties, the queries you ran, and anything about the object
or its surroundings that a careful astronomer would note, including
"nothing notable". Do not assume any of these objects is unusual.
Return one JSON record per position following review_record.schema.json.
