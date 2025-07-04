%This file is part of Directional, a library for directional field processing.
%Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
%
%This Source Code Form is subject to the terms of the Mozilla Public License
%v. 2.0. If a copy of the MPL was not distributed with this file, You can
%obtain one at http://mozilla.org/MPL/2.0/.

%Integration by rounding seams directly.
centers = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
nf=length(F);

rawField(:,4:6)=-rawField(:,4:6);
rawField(:,10:12)=-rawField(:,10:12);

%figure
%hold on
%patch('vertices', V, 'faces', F, 'faceColor', 'w', 'edgeColor', 'b');

%normalizing field to have average length 1 and putting away extra in param length.
avgGradNorm=0;
for i=1:nf
    for j=1:N
        avgGradNorm=avgGradNorm+normv(rawField(i,3*(j-1)+1:3*(j-1)+3));
    end
end
avgGradNorm=avgGradNorm/(N*nf);

rawField=rawField/avgGradNorm;


normals = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:),2);
normals=normals./normv(normals);
%quiver3(centers(:,1), centers(:,2), centers(:,3), normals(:, 1),normals(:, 2),normals(:, 3));

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


rawField = G*x0;
rawField=reshape(rawField, 3*N, length(F))';
funcValues = x2CornerMat*x0;
funcValues = reshape(funcValues, N, length(V))';

%for i=0:N-1
%    axis equal
%    cameratoolbar
%    quiver3(centers(:,1), centers(:,2), centers(:,3), rawField(:, i*3+1),rawField(:, i*3+2),rawField(:, i*3+3));
%end

%testing G2
rawField2 = G2*x0;
rawField2=reshape(rawField2, 2*N, length(F))';
for i=0:N-1
    rawFieldFrom2(:,i*3+1:i*3+3) = rawField2(:,2*i+1).*B1 +rawField2(:,2*i+2).*B2;
end
G2Error = max(max(abs(rawFieldFrom2-rawField)))

funcValues = x2CornerMat*x0;
funcValues = reshape(funcValues, N, length(V))';

%for i=0:N-1
%    axis equal
%    cameratoolbar
%    quiver3(centers(:,1), centers(:,2), centers(:,3), rawField(:, i*3+1),rawField(:, i*3+2),rawField(:, i*3+3));
%end

%figure
%hold on
%patch('vertices', funcValues(:,1:2), 'faces', F, 'faceColor', 'w', 'edgeColor','b');

fixedIndices=[];
fixedValues=[];

leftIndices=integerIndices;
wobj=1;
wconst=10e5;
wbarrier=0.0001;
wclose = 0.01;
s=0.1;
xcurr=x0;
fraction=1;

IImagField=repmat((1:N*nf)',1,4);
JImagField=IImagField;
for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    JImagField(barOffset:barOffset+N-1,1:2)=reshape(varOffset:varOffset+2*N-1, 2, N)';
    JImagField(barOffset:barOffset+N-1,3:4)=JImagField([barOffset+1:barOffset+N-1,barOffset],1:2);
end

origFieldVolumes = zeros(nf,N);
field = G2*x0;
for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(field(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    %tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    origFieldVolumes(i+1,:) = abs((faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1)));
end
i=1;
packetSize=1
roundCoeff=1
success=1;
tic
while (i<=length(integerIndices))
    i
    xprev=xcurr;
    roundDiffs= abs(fraction*xcurr(leftIndices)-roundCoeff*round(fraction*xcurr(leftIndices)/roundCoeff));
    [~,minRoundIndices]=sort(roundDiffs);
    topIndex=min([packetSize;length(minRoundIndices)]);
    minRoundIndices=minRoundIndices(1:topIndex);
    origValues = xcurr(leftIndices(minRoundIndices));
    roundValues = roundCoeff*round(fraction*xcurr(leftIndices(minRoundIndices))/roundCoeff)/fraction;
    prevFixedIndices = fixedIndices;
    prevFixedValues=fixedValues;
    prevLeftIndices=leftIndices;
    fixedIndices=[fixedIndices;leftIndices(minRoundIndices)];
    fixedValues=[fixedValues;roundValues];
    leftIndices(minRoundIndices)=[];
    roundDiffs= abs(fraction*xcurr(fixedIndices)-roundCoeff*round(fraction*xcurr(fixedIndices)/roundCoeff));
    if (max(roundDiffs)<10e-7) %no extra rounding is needed
        i=i+packetSize;
        continue
    end
    options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','display', 'iter-detailed','MaxIterations', 250, 'SpecifyObjectiveGradient',true,'FiniteDifferenceType', 'central');%,'CheckGradients', true); %
    objfun = @(xcurr)objective(xcurr,xprev,A,b,C, G2, N, FN, fixedIndices, fixedValues, s, wobj,wclose, wconst, wbarrier, IImagField, JImagField,origFieldVolumes);
    
    %for checking gradients only!!!
    %x0=rand(size(x0));
    [xcurr ,~,~,exitflag,~]= lsqnonlin(objfun,xcurr,[],[],options);
    roundDiffs= abs(fraction*xcurr(fixedIndices)-roundCoeff*round(fraction*xcurr(fixedIndices)/roundCoeff));
    if (max(roundDiffs)>10e-7) || (size(C,1)~=0 && max(abs(C*xcurr))>10e-7)  %did not converge, try again with different integers
        if (packetSize>1)
            packetSize=floor(packetSize/2)
        else
            fraction=fraction*2
            if (fraction==8)
                success=-1;
                break;
            end
            packetSize=1
        end
        fixedIndices=prevFixedIndices;
        fixedValues=prevFixedValues;
        leftIndices=prevLeftIndices;
        xcurr=xprev;
    else
        i=i+packetSize;
    end
end
toc

xcurr=xcurr*fraction; %so that we do get full integers at the end
xcurr(integerIndices)=round(xcurr(integerIndices));

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

%figure
%hold on
%patch('vertices', funcValues(:,1:2), 'faces', F, 'faceColor', 'none', 'edgeColor','b');
%patch('vertices', funcValuesCurr(:,1:2), 'faces', F, 'faceColor', 'none', 'edgeColor','r');

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



function [f,g]=objective(xcurr,xprev,A,b,C,G2,  N, FN, fixedIndices, fixedValues, s, wobj,wclose, wconst, wbarrier, IImagField, JImagField, origFieldVolumes)

fAb = A*xcurr-b;
gAb = A;
fClose = xcurr-xprev;
gClose=speye(length(xcurr));
fLinConst=C*xcurr;
fConst = xcurr(fixedIndices)-fixedValues;

nf = length(FN);
field = G2*xcurr;
fBarrier = zeros(N*nf,1);
fBarrierDerivative= zeros(N*nf,1);

SImagField=IImagField;

for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(field(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    %tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    imagProduct = (faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1))./origFieldVolumes(i+1,:)';
    barResult = (imagProduct/s).^3 - 3*(imagProduct/s).^2 + 3*(imagProduct/s);
    barResult2 = 1./barResult -1;
    barResult2(imagProduct<=0)=Inf;
    barResult2(imagProduct>=s)=0;
    fBarrier(barOffset:barOffset+N-1)=barResult2;
    
    barDerivative=(3*(imagProduct.^2/s^3) -6*(imagProduct/s^2) + 3/s)./origFieldVolumes(i+1,:)';
    barDerivative(imagProduct<=0)=Inf;
    barDerivative(imagProduct>=s)=0;
    fBarrierDerivative(barOffset:barOffset+N-1) =  barDerivative;
    
    %ImagField(barOffset:barOffset+N-1,1:2)=reshape(varOffset:varOffset+2*N-1, 2, N)';
    %JImagField(barOffset:barOffset+N-1,3:4)=JImagField([barOffset+1:barOffset+N-1,barOffset],1:2);
    
    SImagField(barOffset:barOffset+N-1,:)=[faceFieldNext(:, 2), -faceFieldNext(:, 1), -faceField(:,2),faceField(:,1)];
end

if (nargout<2) %don't compute jacobian
    f=[wobj*fAb;wclose*fClose;wconst*fConst;wconst*fLinConst;wbarrier*fBarrier];
    return
end
gLinConst=C;
gConst=sparse((1:length(fixedIndices))', fixedIndices, ones(length(fixedIndices),1), length(fixedIndices),length(xcurr));
gImagField=sparse(IImagField, JImagField, SImagField, N*nf, 2*N*nf);
barDerVec=-fBarrierDerivative./(fBarrier.^2);
barDerVec(isinf(barDerVec))=0;
barDerVec(isnan(barDerVec))=0;
gBarrierFunc = spdiags(barDerVec, 0, length(fBarrier), length(fBarrier));   %./fBarrier.^2
gBarrier=gBarrierFunc*gImagField*G2;
%f=sum(wobj*fAb.^2)+sum(wconst*fConst.^2)+sum(wbarrier*fBarrier.^2);
%g=2*wobj*gAb'*fAb+2*wconst*gConst'*fConst+2*wbarrier*gBarrier'*fBarrier;
f=[wobj*fAb;wclose*fClose;wconst*fConst;wconst*fLinConst;wbarrier*fBarrier];
g=[wobj*gAb;wclose*gClose;wconst*gConst;wconst*gLinConst;wbarrier*gBarrier];
end





