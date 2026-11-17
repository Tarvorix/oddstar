# Round 2: Adversarial Re-Verification Prompt

A measurement has been reported for Gaia DR3 <id>: <one-sentence
numerical claim, e.g. "W4 flux is N times the photospheric prediction">.
Treat the claim as likely wrong. Starting from raw archive data, not
from any provided file, try to reproduce the number. Then list every
mundane way the number could arise without anything unusual about the
star, and test each one you can. Report what you reproduced, what you
could not, what you ruled out with evidence, and what remains untested.
Return a JSON record following review_record.schema.json.
