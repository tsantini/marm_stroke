function y = kmaskFunc(x,transp,kmask,dim)

switch transp
    case 'notransp'
        y = kmask.*fftnc(x,dim,0,1);
    case 'transp'
        y = ifftnc(kmask.*x,dim,0,1);
    otherwise
        error('unknown input')

end