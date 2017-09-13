# ElecVesicle
vesicle EHD simulation

* Routine in `Main.m`
* Initial configuration defined in `Submit.m`

### Important notes

In the function `Main.m`, the operator `y=LinOP(x)` has a catch!

1. The unknown `x` consists of data lying on vesicle nodes `X(1:end)`, where the first & last nodes are the same (duplicated).
2. The inextensibility condition is built with data lying on `X(1:end-1)`
3. The Stoke BIE operator is built with data lying on `X(2:end)`
