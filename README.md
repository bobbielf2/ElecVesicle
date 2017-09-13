# ElecVesicle
vesicle EHD simulation

* Routine in `Main.m`
* Initial configuration defined in `Submit.m`

### Important notes

In the function `Main.m`, the operator `y=LinOP(x)` has a catch!

1. Each vesicle is defined as `X(1:end)`, where the first & last points are the same
2. The inextensibility condition is built using vesicle nodes `X(1:end-1)`
3. The Stoke BIE operator is built using vesicle nodes `X(2:end)`
