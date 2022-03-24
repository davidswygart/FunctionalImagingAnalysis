function [dfS,dfB] = makeNull(sPre,sStim,bPre,bStim,target)
% throw out any dangling epochs
ne = min([size(sPre,3),size(bPre,3)]);
sPre = sPre(:,:,1:ne);
sStim = sStim(:,:,1:ne);
bPre = bPre(:,:,1:ne);
bStim = bStim(:,:,1:ne);

% unit 1
s = sStim - sPre;
s1 = s(:,:,1:2:end);
s2 = s(:,:,2:2:end);

% unit 0
n = sPre - bPre;
n1 = n(:,:,1:2:end);
n2 = n(:,:,2:2:end);

% big response
b = bStim - bPre;
b1 = b(:,:,1:2:end);
b2 = b(:,:,2:2:end);

switch target
    case 'real'
        dfB = b1;
        dfS = s1;
    case 'inv'
        dfB = s1;
        dfS = b1;
    case '1/1'
        dfB = s2;
        dfS = s1;
    case '0/1'
        dfB = n1;
        dfS = s1;
    case '1/0'
        dfB = s1;
        dfS = n1;
    case '0/0'
        dfB = n2;
        dfS = n1;
    case '-1/1'
        dfB = s2*-1;
        dfS = s1;
end

% throw out any dangling epochs
ne = min([size(dfB,3),size(dfS,3)]);
dfB = dfB(:,:,1:ne);
dfS = dfS(:,:,1:ne);
end