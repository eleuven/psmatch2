# psmatch2
Mahalanobis and Propensity score Matching

psmatch2 is a <a href="www.stata.com">Stata<a/> module that implements full Mahalanobis matching and a variety of propensity score matching methods to adjust for pre-treatment observable differences between a group of treated and a group of untreated.

It provides three commands: -psmatch2- perform the matching, -pstest- reports balancing, and -psgraph- display support for the propensity score.

## Installation from GitHub

Development version:

```
net install psmatch2, from(https://raw.githubusercontent.com/eleuven/psmatch2/master) replace
```

Released versions:

```
net install psmatch2, from(https://raw.githubusercontent.com/eleuven/psmatch2/master/versions/v4.0.13) replace
net install psmatch2, from(https://raw.githubusercontent.com/eleuven/psmatch2/master/versions/v4.0.12) replace
```

Future releases should add a frozen copy under `versions/vX.Y.Z/` so older versions remain installable with `net install`.
