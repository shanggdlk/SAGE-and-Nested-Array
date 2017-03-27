function AB = kr(A,B)
%KR Khatri-Rao product


[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in kr.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab= A(:,f)*transpose(B(:,f));
   %ab = B(:,f)*A(:,f)';
   %AB(:,f)=ab(:);
   AB(:,f)=reshape(transpose(ab),[],1);
end