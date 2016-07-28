function[w13,Pci13]=y13_0618(deta,pp1,pp3)
a(1)=0;b(1)=1;
l(1)=a(1)+0.382*(b(1)-a(1));
u(1)=a(1)+0.618*(b(1)-a(1));
for k=1:500             %k的取值范围尽量大一些
    p0l(:,:,k)=inv(l(k)*inv(pp1)+(1-l(k))*inv(pp3));
    p0u(:,:,k)=inv(u(k)*inv(pp1)+(1-u(k))*inv(pp3));
    if trace(p0l(:,:,k))>trace(p0u(:,:,k))
       if b(k)-l(k)<=deta
           w13=u(k);break;  
       else
           a(k+1)=l(k);b(k+1)=b(k);l(k+1)=u(k);u(k+1)=a(k+1)+0.618*(b(k+1)-a(k+1));
       end
    else
       if u(k)-a(k)<=deta
          w13=l(k);break;
       else
            a(k+1)=a(k);b(k+1)=u(k);u(k+1)=l(k);l(k+1)=a(k+1)+0.382*(b(k+1)-a(k+1));
       end
    end
end
Pci13=inv(w13*inv(pp1)+(1-w13)*inv(pp3));
