RUFUS calls from the child against the two pedigree-defined parents. All records
from the caller's final VCF are retained; these outputs are not restricted to
structural variants and do not add an independent de novo filter.

The VCF sample labels are mapped from the staged BAM names to the pedigree sample
IDs without changing genotype columns or variant records. The original caller
VCF, command, sample mapping and log remain in a unique persistent work directory
for each invocation. Earlier attempts are never removed automatically.

The summary counts VCF records, including each breakend record separately. ALT
classes are exclusive per record: SNV, MNV, sequence INDEL (any length), SYMBOLIC,
BND, MIXED (different ALT classes), or OTHER. INFO/SVTYPE counts retain the caller's
labels separately. These are record counts, not counts of independent events or
validated de novo structural variants.
