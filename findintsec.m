function [PDctl] = findintsec(PD,ref)





PDctl(:,1) = PD + ref;
PDctl(:,2) = PD - ref;

PDctl(PDctl(:,1) >= 360,1) = PDctl(PDctl(:,1) >= 360,1) - 360;
PDctl(PDctl(:,2) >= 360,2) = PDctl(PDctl(:,2) >= 360,2) - 360;

PDctl(PDctl(:,1) < 0,1) = PDctl(PDctl(:,1) < 0,1) + 360;
PDctl(PDctl(:,2) < 0,2) = PDctl(PDctl(:,2) < 0,2) + 360;
