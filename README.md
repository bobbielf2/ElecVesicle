# ElecVesicle
vesicle EHD simulation

* Routine in `Main.m`
* Initial configuration defined in `Submit.m`

### Important notes

In the function Main.m, the operator y=LinOP(x) has a catch!
1. Each vesicle is defined as X(1:end), where the first & last points are the same
2. Each component in the input x correspond to a vesicle at positions X(1:end-1)
3. Each component in the ouput y correspond to a vesicle at positions X(2:end)
