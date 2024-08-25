function pdf=fxy(xi,yi,xx,yy,dx,dy,N)
    pdf=0.*xx;
    xmin=min(min(xx));xmax=max(max(xx));
    ymin=min(min(yy));ymax=max(max(yy));
    for k=1:N
        if((xi(k)>=xmin&&xi(k)<=xmax)&&(yi(k)>=ymin&&yi(k)<=ymax))
            ii=floor((xi(k)-xmin)/dx)+1;
            jj=floor((yi(k)-ymin)/dy)+1;
            pdf(ii,jj)=pdf(ii,jj)+1;
        end
    end
%     pdf=pdf/N;
end