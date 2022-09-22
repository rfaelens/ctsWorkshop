library(RxODE)

tvSimulationModel <- "
KA=0.126;
Vc=2700;
CL=34.1;
Vp=774;
Q=0.688;
TLag = 0.527;

d/dt(A0) = -KA*A0;
alag(A0) = TLag;
d/dt(A1) = KA*A0 - CL/Vc * A1 - Q/Vc * A1 + Q/Vp*A2;
d/dt(A2) = Q/Vc*A1 - Q/Vp*A2;
CONC = A1 / Vc * 1000;
d/dt(CAUC) = CONC;
"
iivModel <- "
KA=0.126*exp(EKA);
Vc=2700*exp(EVc);
CL=34.1*exp(ECL);
Vp=774;
Q=0.688;
TLag = 0.527;

d/dt(A0) = -KA*A0;
alag(A0) = TLag;
d/dt(A1) = KA*A0 - CL/Vc * A1 - Q/Vc * A1 + Q/Vp*A2;
d/dt(A2) = Q/Vc*A1 - Q/Vp*A2;
CONC = A1 / Vc * 1000;
d/dt(CAUC) = CONC;
"
iivModelMod <- "
KA=0.126*exp(EKA);
Vc=2700*exp(EVc);
CL=34.1*exp(ECL) * MOD;
Vp=774;
Q=0.688;
TLag = 0.527;

d/dt(A0) = -KA*A0;
alag(A0) = TLag;
d/dt(A1) = KA*A0 - CL/Vc * A1 - Q/Vc * A1 + Q/Vp*A2;
d/dt(A2) = Q/Vc*A1 - Q/Vp*A2;
CONC = A1 / Vc * 1000;
d/dt(CAUC) = CONC;
"

buildModel <- function(text) {
  RxODE::RxODE(text)
}


simulate <- function(model, dose, nbr.doses, dosing.interval, observationTimes, omega=NULL, nSub=1L, params=NULL) {
  events <- eventTable() %>%
    add.dosing(dose=dose, nbr.doses=nbr.doses, dosing.interval=dosing.interval) %>%
    add.sampling(observationTimes)
  result <- rxSolve(model, events=events, omega = omega, nSub=nSub, params=params)
  result
}

