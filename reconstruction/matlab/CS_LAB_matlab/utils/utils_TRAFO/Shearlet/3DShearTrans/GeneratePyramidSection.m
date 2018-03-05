function [P,F]=GeneratePyramidSection(cubeSide)
P=cell(1,4);
[P{1}.X,P{1}.Y,P{1}.Z]= GenerateXYZ(cubeSide);
F{1}=PolFreq(cubeSide,P{1}.X,P{1}.Y,P{1}.Z);

P{2}.Y=P{1}.X;
P{2}.X=P{1}.Y;
P{2}.Z=P{1}.Z;

F{2}=PolFreq(cubeSide,P{2}.X,P{2}.Y,P{2}.Z);
P{3}.Y=P{1}.Y;
P{3}.X=P{1}.Z;
P{3}.Z=P{1}.X;
F{3}=PolFreq(cubeSide,P{3}.X,P{3}.Y,P{3}.Z);
