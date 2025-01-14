function Index=DetectImag(x)
Index=sum(sum(abs(imag(x))>0)');
if Index~=0
    Index=1;
end
