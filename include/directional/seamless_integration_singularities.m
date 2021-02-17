%This file is part of Directional, a library for directional field processing.
%Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
%
%This Source Code Form is subject to the terms of the Mozilla Public License
%v. 2.0. If a copy of the MPL was not distributed with this file, You can
%obtain one at http://mozilla.org/MPL/2.0/.


%integration by rounding singularity values and then remainder seames if any
clear x0
load /Users/amirvaxman/PatternsParam/build/Release/poisson.mat 
centers = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
nf=length(F);

paramLength = normv(max(V,[],1)-min(V,[],1))*lengthRatio;

%normalizing field to have average length 1 and putting away extra in param length.
avgGradNorm=0;
for i=1:nf
    for j=1:N
        avgGradNorm=avgGradNorm+normv(rawField(i,3*(j-1)+1:3*(j-1)+3));
    end
end
avgGradNorm=avgGradNorm/(N*nf);

rawField=rawField/avgGradNorm;
paramLength=paramLength/avgGradNorm;


%normals
normals = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:),2);
triAreas = normv(normals)/2;
normals=normals./normv(normals);


%face mass matrices
I3=reshape((1:3*N*nf)', 3*N, nf)';
J3=I3;
Mx3 = sparse(I3, J3, repmat(triAreas, 1, 3*N), 3*N*nf, 3*N*nf);
I2=reshape((1:2*N*nf)', 2*N, nf)';
J2=I2;
Mx2 = sparse(I2, J2, repmat(triAreas, 1, 2*N), 2*N*nf, 2*N*nf);

%face bases
B1 = V(F(:,2),:)-V(F(:,1),:);
B1=B1./normv(B1);
B2 = cross(normals, B1);
B2=B2./normv(B2);



%creating 3D->2D gradient
IReduc = repmat((1:2*N*nf)',1,3);
JReduc = repmat(reshape((1:3*N*nf)',3,N*nf)', 1,2);
JReduc = reshape(JReduc', 3, 2*N*nf)';
B1N = reshape(repmat(B1', N, 1), 3, N*nf);
B2N = reshape(repmat(B2', N, 1), 3, N*nf);
SReduc = reshape([B1N;B2N], 3, 2*N*nf)';

reducMat = sparse(IReduc, JReduc, SReduc, 2*N*nf, 3*N*nf);
G2 = reducMat*G;

%Reducing C into U
[I,J,S]=find(C);
small2Big=unique(J);
nonPartIndices = setdiff((1:size(C,2))', small2Big);
big2Small=zeros(size(C,2),1);
big2Small(small2Big)=1:length(small2Big);

CSmall = zeros(size(C,1), size(small2Big,1));
CSmall(sub2ind(size(CSmall), I, big2Small(J)))=S;

USmall=null(CSmall);

URaw = blkdiag(speye(length(nonPartIndices)), USmall);
%parmutation matrix to build full uncompressed variables
permMat = sparse([nonPartIndices;small2Big], (1:size(URaw,1))', ones(size(URaw,1),1), size(C,2), size(URaw,1));

UFull = permMat*URaw;

%figure
%hold on
%axis equal
%cameratoolbar
%for i=0:N-1
%    axis equal
%    cameratoolbar
%    quiver3(centers(:,1), centers(:,2), centers(:,3), rawField(:, i*3+1),rawField(:, i*3+2),rawField(:, i*3+3));
%end

%gradientField = G*x0;
%gradientField=reshape(gradientField, 3*N, length(F))';
%figure
%hold on
%axis equal
%cameratoolbar
%for i=0:N-1
%    axis equal
%    cameratoolbar
%    quiver3(centers(:,1), centers(:,2), centers(:,3), gradientField(:, i*3+1),gradientField(:, i*3+2),gradientField(:, i*3+3));
%end


%rawField=reshape(rawField, 3*N, length(F))';
%funcValues = x2CornerMat*x0;
%funcValues = reshape(funcValues, N, length(V))';

%for i=0:N-1
%    axis equal
%    cameratoolbar
%    quiver3(centers(:,1), centers(:,2), centers(:,3), rawField(:, i*3+1),rawField(:, i*3+2),rawField(:, i*3+3));
%end

%testing G2
%rawField2 = G2*x0;
%initField=rawField2;
%rawField2=reshape(rawField2, 2*N, length(F))';
%for i=0:N-1
%    rawFieldFrom2(:,i*3+1:i*3+3) = rawField2(:,2*i+1).*B1 +rawField2(:,2*i+2).*B2;
%end
%G2Error = max(max(abs(rawFieldFrom2-rawField)))

%funcValues = x2CornerMat*x0;
%funcValues = reshape(funcValues, N, length(V))';

%for i=0:N-1
%    axis equal
%    cameratoolbar
%    quiver3(centers(:,1), centers(:,2), centers(:,3), rawField(:, i*3+1),rawField(:, i*3+2),rawField(:, i*3+3));
%end

%figure
%hold on
%patch('vertices', funcValues(:,1:2), 'faces', F, 'faceColor', 'w', 'edgeColor','b');

%rawField=reshape(rawField, 3*N, length(F))';

rawField2=zeros(length(F), 2*N);
for i=0:N-1
    currField = rawField(:,i*3+1:i*3+3);
    rawField2(:, i*2+1:i*2+2) = [dot(currField, B1,2), dot(currField, B2,2)];
end
%HACKING THE FIELD
%rawComplex = complex(rawField2(:,1), rawField2(:,2))*exp(1i*pi/10);
%rawField2(:,3:4) = [real(rawComplex), imag(rawComplex)];
%rawField2(:,7:8) = -rawField2(:,3:4);

rawField2 = reshape(rawField2', 2*N*length(F),1);
%rawField = reshape(rawField', 3*N*length(F),1);

IImagField=repmat((1:N*nf)',1,4);
JImagField=IImagField;
for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    JImagField(barOffset:barOffset+N-1,1:2)=reshape(varOffset:varOffset+2*N-1, 2, N)';
    JImagField(barOffset:barOffset+N-1,3:4)=JImagField([barOffset+1:barOffset+N-1,barOffset],1:2);
end

origFieldVolumes = zeros(nf,N);

for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(rawField2(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    %tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    crossProducts = (faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1));
    origFieldVolumes(i+1,:) = normv(faceField).*normv(faceFieldNext);%crossProducts';
end

minOrigFieldVolumes=min(origFieldVolumes)


wintegration=10e3;
wconst=10e3;
wbarrier=0.0001;
wclose=1;
s=0.1;
success=1;

%Use if Mx doesn't exist
Mx = speye(2*N*size(F,1));


%%%%%%%%%%%%%%%%%%%%%%solving initial poisson problem without integers

L = G2'*Mx*G2;
E = UFull'*G2'*Mx*G2*UFull;
f = UFull'*G2'*Mx*(rawField2/paramLength);
constMat = UFull(fixedIndices,:);
bigMat = [E constMat'; constMat, sparse([],[],[],size(constMat,1), size(constMat,1))];
bigRhs = [f;fixedValues];
initXSmall = bigMat\bigRhs;
initXSmall=initXSmall(1:size(UFull,2));

%initXSmall = (UFull'*G'*Mx*G*UFull)\(UFull'*G'*Mx*(rawField/paramLength));
initialIntegrationError = max(abs(rawField2 - paramLength*G2*UFull*initXSmall))
x0=UFull*initXSmall;
%initXSmall=UFull\x0;  %should get rid of x0...
%for checking gradients only!!!
%initXSmall = rand(size(initXSmall));

sizeX=length(x0);
%just for test!!!
%rawField2 = G2*UFull*initXSmall;
initXandFieldSmall = [initXSmall;rawField2];
currXandFieldSmall=initXandFieldSmall;
UExt = blkdiag(UFull, speye(length(rawField2)));
%[I,J,S]=find(C);
%CExt=sparse(I,J,S,size(C,1), length(initXandField));
options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','display', 'iter-detailed','MaxIterations', 250, 'SpecifyObjectiveGradient',true,'FiniteDifferenceType', 'central');%,'CheckGradients', true);
objfun = @(currXandFieldSmall)initialPoissonObjective(currXandFieldSmall,rawField2,sizeX,UExt, G2,  N, FN, fixedIndices, fixedValues, paramLength, s, wintegration,wclose, wconst, wbarrier, IImagField, JImagField,origFieldVolumes);

%for checking gradients only!!!
%currXandField(1:sizeX) = rand(sizeX,1);
[currXandFieldSmall ,~,~,exitflag,~]= lsqnonlin(objfun,currXandFieldSmall,[],[],options);

currXandField=UExt*currXandFieldSmall;

x0=currXandField(1:sizeX);
x0Small=currXandFieldSmall(1:length(initXSmall));
rawField2=currXandField(sizeX+1:end);

rawField2=reshape(rawField2, 2*N,nf)';
avgGradNorm=0;
for i=1:nf
    for j=1:N
        avgGradNorm=avgGradNorm+normv(rawField2(i,2*(j-1)+1:2*(j-1)+2));
    end
end
avgGradNorm=avgGradNorm/(N*nf);

rawField2=reshape(rawField2', 2*N*nf,1);
rawField=rawField/avgGradNorm;
rawField2=rawField2/avgGradNorm;
x0=x0/avgGradNorm;

if (size(C,1)~=0 && max(abs(C*x0))>10e-7)  %did not converge
    success=-1;
end
%%%%%%%%%%%%%%%%%%%%%solving integer problem

wpoisson=1;
wclose = 0.01;
wconst=10e5;

origFieldVolumes = zeros(nf,N);
field = rawField2;
for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(field(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    %tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    origFieldVolumes(i+1,:) =normv(faceField).*normv(faceFieldNext);% (faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1));
end


minOrigFieldVolumes=min(origFieldVolumes)

fixedIndices=[];
fixedValues=[];


%should transform integer indices...


leftIndices=singularIndices;
xcurrSmall=x0Small;
fraction=1;
roundingOffset=0;

i=1;
roundCoeff=1


tic
while (i<=length(singularIndices))
    i
    if (success==-1)
        break;
    end
    xprevSmall=xcurrSmall;
    xcurr=UFull*xcurrSmall;
    roundDiffs= abs(fraction*xcurr(leftIndices)-roundingOffset-round(fraction*xcurr(leftIndices)-roundingOffset));
    [~,minRoundIndices]=sort(roundDiffs);
    %topIndex=min([packetSize;length(minRoundIndices)]);
    minRoundIndices=minRoundIndices(1);
    origValues = xcurr(leftIndices(minRoundIndices));
    roundValues = (round(fraction*xcurr(leftIndices(minRoundIndices))-roundingOffset)+roundingOffset)/fraction;
    prevFixedIndices = fixedIndices;
    prevFixedValues=fixedValues;
    prevLeftIndices=leftIndices;
    fixedIndices=[fixedIndices;leftIndices(minRoundIndices)];
    fixedValues=[fixedValues;roundValues];
    leftIndices(minRoundIndices)=[];
    roundDiffs= abs(fraction*xcurr(fixedIndices)-roundingOffset-round(fraction*xcurr(fixedIndices)-roundingOffset));
    if (max(roundDiffs)<10e-7) %no extra rounding is needed
        i=i+1;
        continue
    end
    options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','display', 'iter-detailed','MaxIterations', 250, 'SpecifyObjectiveGradient',true,'FiniteDifferenceType', 'central');%,'CheckGradients', true); %
    objfun = @(xcurrSmall)integerPoissonObjective(xcurrSmall,xprevSmall,rawField2, UFull, G2, Mx, N, FN, fixedIndices, fixedValues, paramLength, s, wpoisson,wclose, wconst, wbarrier, IImagField, JImagField,origFieldVolumes);
    
    %for checking gradients only!!!
    %x0=rand(size(x0));
    [xcurrSmall ,~,~,exitflag,~]= lsqnonlin(objfun,xcurrSmall,[],[],options);
    xcurr = UFull*xcurrSmall;
    roundDiffs= abs(fraction*xcurr(fixedIndices)-roundingOffset-round(fraction*xcurr(fixedIndices)-roundingOffset));
    if (max(roundDiffs)>10e-7) || (size(C,1)~=0 && max(abs(C*xcurr))>10e-7)  %did not converge, try again with different integers
        %fraction=fraction*2
        %if (fraction==8)
            success=-1;
            break;
        %end
        %fixedIndices=prevFixedIndices;
        %fixedValues=prevFixedValues;
        %leftIndices=prevLeftIndices;
        %xcurrSmall=xprevSmall;
    else
        i=i+1;
    end
end
toc

xcurr = UFull*xcurrSmall;
xcurr=xcurr*fraction; %so that we do get full integers at the end
xcurr(singularIndices)=round((xcurr(singularIndices)));
integerError = max(abs(xcurr(integerIndices)-round(xcurr(integerIndices))))
if (integerError>10e-7)
    tic
    i=1;
    fixedIndices=[];
    fixedValues=[];
    leftIndices=integerIndices;
    while (i<=length(integerIndices))
        i
        if (success==-1)
            break;
        end
        xprevSmall=xcurrSmall;
        xcurr=UFull*xcurrSmall;
        roundDiffs= abs(fraction*xcurr(leftIndices)-roundingOffset-round(fraction*xcurr(leftIndices)-roundingOffset));
        [~,minRoundIndices]=sort(roundDiffs);
        %topIndex=min([packetSize;length(minRoundIndices)]);
        minRoundIndices=minRoundIndices(1);
        origValues = xcurr(leftIndices(minRoundIndices));
        roundValues = (round(fraction*xcurr(leftIndices(minRoundIndices))-roundingOffset)+roundingOffset)/fraction;
        prevFixedIndices = fixedIndices;
        prevFixedValues=fixedValues;
        prevLeftIndices=leftIndices;
        fixedIndices=[fixedIndices;leftIndices(minRoundIndices)];
        fixedValues=[fixedValues;roundValues];
        leftIndices(minRoundIndices)=[];
        roundDiffs= abs(fraction*xcurr(fixedIndices)-roundingOffset-round(fraction*xcurr(fixedIndices)-roundingOffset));
        if (max(roundDiffs)<10e-7) %no extra rounding is needed
            i=i+1;
            continue
        end
        options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','display', 'iter-detailed','MaxIterations', 250, 'SpecifyObjectiveGradient',true,'FiniteDifferenceType', 'central');%,'CheckGradients', true); %
        objfun = @(xcurrSmall)integerPoissonObjective(xcurrSmall,xprevSmall,rawField2, UFull, G2, Mx, N, FN, fixedIndices, fixedValues, paramLength, s, wpoisson,wclose, wconst, wbarrier, IImagField, JImagField,origFieldVolumes);
        
        %for checking gradients only!!!
        %x0=rand(size(x0));
        [xcurrSmall ,~,~,exitflag,~]= lsqnonlin(objfun,xcurrSmall,[],[],options);
        xcurr = UFull*xcurrSmall;
        roundDiffs= abs(fraction*xcurr(fixedIndices)-roundingOffset-round(fraction*xcurr(fixedIndices)-roundingOffset));
        if (max(roundDiffs)>10e-7) || (size(C,1)~=0 && max(abs(C*xcurr))>10e-7)  %did not converge, try again with different integers
            %fraction=fraction*2
            %if (fraction==8)
                success=-1;
                break;
            %end
            %fixedIndices=prevFixedIndices;
            %fixedValues=prevFixedValues;
            %leftIndices=prevLeftIndices;
            %xcurrSmall=xprevSmall;
        else
            i=i+1;
        end
    end
    toc
end

xcurr = UFull*xcurrSmall;
xcurr=xcurr*fraction; %so that we do get full integers at the end
integerError = max(abs(xcurr(integerIndices)-round(xcurr(integerIndices))))
xcurr(integerIndices)=round(xcurr(integerIndices));

%adding roundingoffset
%xcurr(setdiff(1:length(xcurr), integerIndices))=xcurr(setdiff(1:length(xcurr), integerIndices))+0.5;


rawFieldCurr = G*xcurr;
%rawFieldCurr=reshape(rawFieldCurr, 3*N, length(F))';
funcValuesCurr = x2CornerMat*xcurr;
funcValuesCurr = reshape(funcValuesCurr, N, length(V))';

%checking bijectivity
minBijectivity=zeros(length(FN),1);
for i=0:length(FN)-1
    varOffset=i*3*N+1;
    barOffset=i*N+1;
    faceField = reshape(rawFieldCurr(varOffset:varOffset+3*N-1),3,N)';
    faceField=faceField./repmat(normv(faceField), 1,3);
    tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    minBijectivity(i+1)=min(tripleProducts);
end

figure
hold on
patch('vertices', funcValuesCurr(:,1:2), 'faces', F, 'faceColor', 'none', 'edgeColor','r');
axis equal


%figure
%hold on
%patch('vertices', funcValuesCurr(:,1:2), 'faces', F, 'faceColor', 'flat', 'CData', minBijectivity);

%for i=0:N-1
%    figure
%    hold on
%    patch('vertices', V, 'faces', F, 'faceColor', 'interp', 'CData', funcValues(:,i+1));
%   axis equal
%    cameratoolbar
%    quiver3(centers(:,1), centers(:,2), centers(:,3), rawField(:, i*3+1),rawField(:, i*3+2),rawField(:, i*3+3));
%end

function [f,g]=initialPoissonObjective(xAndCurrFieldSmall,initField,sizeX, UExt,G2, N, FN, fixedIndices, fixedValues, paramLength, s, wintegration,wclose, wconst, wbarrier, IImagField, JImagField, origFieldVolumes)

xAndCurrField = UExt*xAndCurrFieldSmall;
xcurr=xAndCurrField(1:sizeX);
currField=xAndCurrField(sizeX+1:end);

fIntegration = (currField - paramLength*G2*xcurr);
gIntegration = [-paramLength*G2, speye(length(currField))]*UExt;

fClose = currField-initField;
gClose = sparse(1:length(currField), sizeX+1:length(xAndCurrField), ones(1,length(currField)), length(currField), length(xAndCurrField))*UExt;

%fLinConst=CExt*xAndCurrField;
%gLinConst=CExt;

fConst = xcurr(fixedIndices)-fixedValues;
gConst = sparse((1:length(fixedIndices))', fixedIndices, ones(length(fixedIndices),1), length(fixedIndices),length(xAndCurrField))*UExt;

nf = length(FN);

fBarrier = zeros(N*nf,1);
splineDerivative= zeros(N*nf,1);
barSpline=zeros(N*nf,1);

SImagField=IImagField;

for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(currField(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    %tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    imagProduct = (faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1))./origFieldVolumes(i+1,:)';
    barResult = (imagProduct/s).^3 - 3*(imagProduct/s).^2 + 3*(imagProduct/s);
    barResult2 = 1./barResult -1;
    
    barResult2(imagProduct<=0)=Inf;
    barResult2(imagProduct>=s)=0;
    fBarrier(barOffset:barOffset+N-1)=barResult2;
    barSpline(barOffset:barOffset+N-1)=barResult;
    
    splineDerivativeLocal=(3*(imagProduct.^2/s^3) -6*(imagProduct/s^2) + 3/s);
    splineDerivativeLocal(imagProduct<=0)=Inf;
    splineDerivativeLocal(imagProduct>=s)=0;
    splineDerivative(barOffset:barOffset+N-1) =  splineDerivativeLocal;
    
    %ImagField(barOffset:barOffset+N-1,1:2)=reshape(varOffset:varOffset+2*N-1, 2, N)';
    %JImagField(barOffset:barOffset+N-1,3:4)=JImagField([barOffset+1:barOffset+N-1,barOffset],1:2);
    
    SImagField(barOffset:barOffset+N-1,:)=[faceFieldNext(:, 2), -faceFieldNext(:, 1), -faceField(:,2),faceField(:,1)]./origFieldVolumes(i+1,:)';
end

if (nargout<2) %don't compute jacobian
    f=[wintegration*fIntegration;wclose*fClose;wconst*fConst;wbarrier*fBarrier];
    %f=fBarrier;
    return
end


gImagField=sparse(IImagField, JImagField, SImagField, N*nf, length(currField));
gFieldReduction = gClose;
barDerVec=-splineDerivative./(barSpline.^2);
barDerVec(fBarrier==Inf)=Inf;
barDerVec(isinf(barDerVec))=0;
barDerVec(isnan(barDerVec))=0;
gBarrierFunc = spdiags(barDerVec, 0, length(fBarrier), length(fBarrier));   %./fBarrier.^2
gBarrier=gBarrierFunc*gImagField*gFieldReduction;
%f=sum(wobj*fAb.^2)+sum(wconst*fConst.^2)+sum(wbarrier*fBarrier.^2);
%g=2*wobj*gAb'*fAb+2*wconst*gConst'*fConst+2*wbarrier*gBarrier'*fBarrier;
f=[wintegration*fIntegration;wclose*fClose;wconst*fConst;wbarrier*fBarrier];
g=[wintegration*gIntegration;wclose*gClose;wconst*gConst;wbarrier*gBarrier];
%f=fBarrier;
%g=gBarrier;
end




function [f,g]=integerPoissonObjective(xcurrSmall,xprevSmall,rawField2, U, G2, sqrtMx,  N, FN, fixedIndices, fixedValues, paramLength, s, wobj,wclose, wconst, wbarrier, IImagField, JImagField, origFieldVolumes)

xcurr = U*xcurrSmall;
fObj = G2*U*xcurrSmall*paramLength - rawField2;
gObj = G2*U*paramLength;
%fObj = A*xcurr-b;
%gObj = A*U;
fClose = (xcurrSmall-xprevSmall);
gClose=speye(length(xcurrSmall));
%fLinConst=C*xcurr;

fConst = (xcurr(fixedIndices)-fixedValues);

nf = length(FN);
currField = G2*xcurr*paramLength;
fBarrier = zeros(N*nf,1);
splineDerivative= zeros(N*nf,1);
barSpline=zeros(N*nf,1);

SImagField=IImagField;

for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(currField(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    %tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    imagProduct = (faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1))./origFieldVolumes(i+1,:)';
    barResult = (imagProduct/s).^3 - 3*(imagProduct/s).^2 + 3*(imagProduct/s);
    barResult2 = 1./barResult -1;
    
    barResult2(imagProduct<=0)=Inf;
    barResult2(imagProduct>=s)=0;
    fBarrier(barOffset:barOffset+N-1)=barResult2;
    barSpline(barOffset:barOffset+N-1)=barResult;
    
    splineDerivativeLocal=(3*(imagProduct.^2/s^3) -6*(imagProduct/s^2) + 3/s);
    splineDerivativeLocal(imagProduct<=0)=Inf;
    splineDerivativeLocal(imagProduct>=s)=0;
    splineDerivative(barOffset:barOffset+N-1) =  splineDerivativeLocal;
    
    %ImagField(barOffset:barOffset+N-1,1:2)=reshape(varOffset:varOffset+2*N-1, 2, N)';
    %JImagField(barOffset:barOffset+N-1,3:4)=JImagField([barOffset+1:barOffset+N-1,barOffset],1:2);
    
    SImagField(barOffset:barOffset+N-1,:)=[faceFieldNext(:, 2), -faceFieldNext(:, 1), -faceField(:,2),faceField(:,1)]./origFieldVolumes(i+1,:)';
end

if (nargout<2) %don't compute jacobian
    f=[wobj*fObj;wclose*fClose;wconst*fConst;wbarrier*fBarrier];
    return
end
%gLinConst=C;
gConst=sparse((1:length(fixedIndices))', fixedIndices, ones(length(fixedIndices),1), length(fixedIndices),length(xcurr))*U;
gImagField=sparse(IImagField, JImagField, SImagField, N*nf, length(currField));
barDerVec=-splineDerivative./(barSpline.^2);
barDerVec(fBarrier==Inf)=Inf;
barDerVec(isinf(barDerVec))=0;
barDerVec(isnan(barDerVec))=0;
gBarrierFunc = spdiags(barDerVec, 0, length(fBarrier), length(fBarrier));   %./fBarrier.^2
gBarrier=gBarrierFunc*gImagField*G2*U*paramLength;
%f=sum(wobj*fAb.^2)+sum(wconst*fConst.^2)+sum(wbarrier*fBarrier.^2);
%g=2*wobj*gAb'*fAb+2*wconst*gConst'*fConst+2*wbarrier*gBarrier'*fBarrier;
f=[wobj*fObj;wclose*fClose;wconst*fConst;wbarrier*fBarrier];
g=[wobj*gObj;wclose*gClose;wconst*gConst;wbarrier*gBarrier];
end





