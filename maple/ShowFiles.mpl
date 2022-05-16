ScanOneEq:=proc(lines)
    # Scans one set of equations (input is the set of all lines between "Subject To" and "END")
    local nx, str, i, tmp, j, res;
    res:=[]:
    for i from 1 to nops(lines) do
        nx:=0: # number of monomials in the equation
        str:="  %*c%*c%*d:";
        tmp:=sscanf(lines[i],cat(str, " %d %s"));
        while ((nops(tmp) = 2*(nx+1)) and (substring(tmp[-1],1) = "x")) do
            str:=cat(str," %d %s");
            nx:=nx+1;
            tmp:=sscanf(lines[i],cat(str, " %d %s"));
        od;
        str:=cat(str," %[>=] %d");
        tmp:=sscanf(lines[i],str);
        if (tmp[-2] = ">=") then
            res:=[op(res),add(tmp[2*j-1]*convert(tmp[2*j],name),j=1..nops(tmp)/2-1)>=tmp[-1]];
        else
            res:=[op(res),add(tmp[2*j-1]*convert(tmp[2*j],name),j=1..nops(tmp)/2-1)=tmp[-1]];
        fi;
    od;
    return res;
end:



ScanFile:=proc(file)
    # Scans all equations provided in a given file
    local line,lines,res,tmp;
    tmp:=1;
    line:=readline(file);
    res:=[];
    while (line <> 0) do
        while (line <> "Subject To") do
            line:=readline(file):
            if (line = 0) then
                return res
            fi
        od:
        lines:=[readline(file)]:
        while ((lines[-1] <> "END") and (lines[-1] <> "BOUNDS")) do
            lines:=[op(lines),readline(file)]:
        od:
        tmp:=tmp+1;
        res:=[op(res),ScanOneEq(lines[1..-2])];
        line:=readline(file);
    od;
    return res;
end:



MakePolyhedra:=proc(file)
    # Make the maple polyhedra of the tropical branches described in the file.
    local res, relation, branch, branchtab, i , varset, varnum, coord;
    with(PolyhedralSets);
    varset:={};
    branchtab:=[];
    #pbranchtab:=[];
    res :=ScanFile(file); # on vient de récupérer le gros tableau de réponses de ptcut en format maple : [[polyt1],...,[polytn]].
    coord:=sort(Coordinates(PolyhedralSet(res[1])),length);
    for i from 1 to numelems(res) do
        relation:=res[i];
        branch:=PolyhedralSet(relation,coord);
        branchtab:=[op(branchtab),branch];
        #pbranch:=Project(branch,[proj1,proj2,proj3],reducespace);
        #pbranchtab:=[op(pbranchtab),pbranch];
    od;
    return branchtab;
end:



Show3DLogTraj:=proc(file,eps,nbrcolumn,listnotread)
    # Show trajectories of the columns to be read in the given datafile in log[eps] basis. nbrcolumn is the number of columns in the given datafile, and listnotread is a list of columns to be deleted from the data. It should remain 3 columns.
    local d, dd, ddd, logtraj, B, A, C, n, i, l;
    with(plots);
    d:=readdata(file,nbrcolumn);
    ddd:=LinearAlgebra[DeleteColumn](d,listnotread);
    #n:=LinearAlgebra[RowDimension](ddd);
    #i:=1;
    #l:=[];
    #while i+9<n do
    #	l:=[op(l),i,i+1,i+3,i+4,i+5,i+6,i+7,i+8,i+9];
    #	i:=i+10;
    #od;
    dd:=LinearAlgebra[DeleteRow](ddd,1);
    logtraj:=map(log[eps],dd);
    B:=plots[textplot3d]([[logtraj[1,1],logtraj[1,2],logtraj[1,3],"alpha",color="Red"],[logtraj[-1,1],logtraj[-1,2],logtraj[-1,3],"omega",color="Blue"]]);
    A:=plots[pointplot3d](logtraj,color="Black",symbol=solidsphere,symbolsize=2);
    C:=plots[display](A,B);
    return C;
end:



Show3DBranches:=proc(file,proj1,proj2,proj3,name1,name2,name3)
	# Show the projection of the tropical branches described in the file on [proj1,proj2,proj3].
    local figure, res, relation, branch, pbranch, branchtab, pbranchtab, i , k, kstring, varset, varnum, colorlist, labels0;
    with(PolyhedralSets);
    colorlist:=[COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63)];
    varset:={};
    branchtab:=[];
    pbranchtab:=[];
    res :=ScanFile(file); # on vient de récupérer le gros tableau de réponses de ptcut en format maple : [[polyt1],...,[polytn]].
    for i from 1 to numelems(res) do
        relation:=res[i];
        branch:=PolyhedralSet(relation);
        branchtab:=[op(branchtab),branch];
        pbranch:=Project(branch,[proj1,proj2,proj3],reducespace);
        pbranchtab:=[op(pbranchtab),pbranch];
    od;
    labels0:=[name1,name2,name3];
    figure:=Plot(pbranchtab, labels=labels0,labeldirections=[horizontal,horizontal,horizontal], color=colorlist[1..numelems(res)], transparency=0.5);
    #with(plots);
    #plots[display](figure);
    return figure;
end:



Show3DBranchesWithTraj:=proc(file1,file2,eps,nbcolumn,listnotread,proj1,proj2,proj3,name1,name2,name3)
	# Show the projection of the tropical branches described in the file1 on [proj1,proj2,proj3] and the log[eps] trajectories of file2.
    local logtrajfig, fig, figure, res, relation, branch, pbranch, branchtab, pbranchtab, i , k, kstring, varset, varnum, colorlist, labels0;
    with(PolyhedralSets);
    with(plots);
    colorlist:=[COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63)];
    varset:={};
    branchtab:=[];
    pbranchtab:=[];
    res :=ScanFile(file1); # on vient de récupérer le gros tableau de réponses de ptcut en format maple : [[polyt1],...,[polytn]].
    for i from 1 to numelems(res) do
        relation:=res[i];
        branch:=PolyhedralSet(relation);
        branchtab:=[op(branchtab),branch];
        pbranch:=Project(branch,[proj1,proj2,proj3],reducespace);
        pbranchtab:=[op(pbranchtab),pbranch];
    od;
    labels0:=[name1,name2,name3];
    fig:=Plot(pbranchtab, labels=labels0,labeldirections=[horizontal,horizontal,horizontal], color=colorlist[1..numelems(res)], transparency=0.5);
    logtrajfig:=Show3DLogTraj(file2,eps,nbcolumn,listnotread);
    figure:=plots[display](fig,logtrajfig);
    return figure;
end:



Show3DBranchesFromList:=proc(list,proj1,proj2,proj3,name1,name2,name3)
	#Show the projection of the tropical branches described in the list (generally obtained with MakePolyhedra) on [proj1,proj2,proj3]
	local figure, res, relation, branch, pbranch, branchtab, pbranchtab, i , k, kstring, varset, varnum, colorlist, labels0;
	with(PolyhedralSets);
    with(plots);
    colorlist:=[COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63)];
    varset:={};
    branchtab:=[];
    pbranchtab:=[];
    for i from 1 to numelems(list) do
    	pbranch:=Project(list[i],[proj1,proj2,proj3],reducespace);
    	pbranchtab:=[op(pbranchtab),pbranch];
    od;
    labels0:=[name1,name2,name3];
    figure:=Plot(pbranchtab, labels=labels0,labeldirections=[horizontal,horizontal,horizontal], color=colorlist[1..numelems(list)], transparency=0.5);
    return figure;
end:



Show3DBranchesFromListWithTraj:=proc(list,file2,eps,nbcolumn,listnotread,proj1,proj2,proj3,name1,name2,name3)
	#Show the projection of the tropical branches described in the list (generally obtained with MakePolyhedra) on [proj1,proj2,proj3] and the log[eps] trajectories of file2.
	local logtrajfig, fig, figure, res, relation, branch, pbranch, branchtab, pbranchtab, i , k, kstring, varset, varnum, colorlist, labels0;
	with(PolyhedralSets);
    with(plots);
    colorlist:=[COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63)];
    varset:={};
    branchtab:=[];
    pbranchtab:=[];
    for i from 1 to numelems(list) do
    	pbranch:=Project(list[i],[proj1,proj2,proj3],reducespace);
    	pbranchtab:=[op(pbranchtab),pbranch];
    od;
    labels0:=[name1,name2,name3];
    fig:=Plot(pbranchtab, labels=labels0,labeldirections=[horizontal,horizontal,horizontal], color=colorlist[1..numelems(list)], transparency=0.5);
    logtrajfig:=Show3DLogTraj(file2,eps,nbcolumn,listnotread);
    figure:=plots[display](fig,logtrajfig);
end:



Show2DBranches:=proc(file,proj1,proj2,name1,name2,transp:=0.5,thickn:=1)
	#Show the projection of the tropical branches described in the file on [proj1,proj2].
    local figure, res, relation, branch, pbranch, branchtab, pbranchtab, i , k, kstring, varset, varnum, colorlist, labels0;
    with(PolyhedralSets);
    #colorlist:=[COLOUR(RGB,255,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,255), COLOUR(RGB,0,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63)];
    colorlist:=[COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63)];
    varset:={};
    branchtab:=[];
    pbranchtab:=[];
    res :=ScanFile(file); # on vient de récupérer le gros tableau de réponses de ptcut en format maple : [[polyt1],...,[polytn]].
    for i from 1 to numelems(res) do
        relation:=res[i];
        branch:=PolyhedralSet(relation);
        branchtab:=[op(branchtab),branch];
        pbranch:=Project(branch,[proj1,proj2],reducespace);
        pbranchtab:=[op(pbranchtab),pbranch];
    od;
    labels0:=[name1,name2];
    figure:=Plot(pbranchtab, labels=labels0,labeldirections=[horizontal,vertical], color=colorlist[1..numelems(res)], thickness=thickn, transparency=transp);#, vertexoptions=[symbolsize=1]);
    #with(plots);
    #plots[display](figure);
    return figure;
end:



Show2DLogTraj:=proc(file,eps,nbrcolumn,listnotread,name1,name2)
    #Show trajectories of the columns to be read in the given datafile in log[eps] basis. nbrcolumn is the number of columns in the given datafile, and listnotread is a list of columns to be deleted from the data. It should remain 2 columns.
    local d, dd, ddd, logtraj, B, A, C, labels0;
    with(plots);
    d:=readdata(file,nbrcolumn);
    ddd:=LinearAlgebra[DeleteColumn](d,listnotread);
    dd:=LinearAlgebra[DeleteRow](ddd,1);
    logtraj:=map(log[eps],dd);
    labels0:=[name1,name2];
    B:=plots[textplot]([[logtraj[1,1],logtraj[1,2],"alpha",color="Red"],[logtraj[-1,1],logtraj[-1,2],"omega",color="Blue"]]);
    A:=plots[pointplot](logtraj,color="Black",symbolsize=5,labels=labels0);
    #A:=plots[pointplot](logtraj,color="Red",symbol=solidsphere,symbolsize=5);
    C:=plots[display](A,B);
    return C;
end:



Show2DBranchesWithTraj:=proc(file1,file2,eps,nbcolumn,listnoread,proj1,proj2,name1,name2)
	#Show the projection of the tropical branches described in the file on [proj1,proj2] with the log[eps] trajectories of file2.
    local logtrajfig, fig, figure, res, relation, branch, pbranch, branchtab, pbranchtab, i , k, kstring, varset, varnum, colorlist, labels0;
    with(PolyhedralSets);
    colorlist:=[COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63), COLOUR(RGB,255,0,0), COLOUR(RGB,0,255,0), COLOUR(RGB,0,0,255), COLOUR(RGB,255,255,0), COLOUR(RGB,255,0,255), COLOUR(RGB,0,255,255), COLOUR(RGB,127,127,127), COLOUR(RGB,63,255,0), COLOUR(RGB,0,63,255), COLOUR(RGB,255,0,63), COLOUR(RGB,0,127,127), COLOUR(RGB,127,0,127), COLOUR(RGB,127,127,0), COLOUR(RGB,0,193,255), COLOUR(RGB,255,0,193), COLOUR(RGB,193,255,0), COLOUR(RGB,0,63,193), COLOUR(RGB,193,0,63), COLOUR(RGB,63,193,0), COLOUR(RGB,255,127,127), COLOUR(RGB,127,255,127), COLOUR(RGB,127,127,255), COLOUR(RGB,193,0,0), COLOUR(RGB,0,193,0), COLOUR(RGB,0,0,193), COLOUR(RGB,255,255,63), COLOUR(RGB,255,63,255), COLOUR(RGB,63,255,255), COLOUR(RGB,193,63,0), COLOUR(RGB,0,193,63), COLOUR(RGB,63,0,193), COLOUR(RGB,63,127,0), COLOUR(RGB,0,63,127), COLOUR(RGB,127,0,63), COLOUR(RGB,63,255,193), COLOUR(RGB,193,63,255), COLOUR(RGB,255,193,63)];
    varset:={};
    branchtab:=[];
    pbranchtab:=[];
    res :=ScanFile(file); # on vient de récupérer le gros tableau de réponses de ptcut en format maple : [[polyt1],...,[polytn]].
    for i from 1 to numelems(res) do
        relation:=res[i];
        branch:=PolyhedralSet(relation);
        branchtab:=[op(branchtab),branch];
        pbranch:=Project(branch,[proj1,proj2],reducespace);
        pbranchtab:=[op(pbranchtab),pbranch];
    od;
    labels0:=[name1,name2];
    fig:=Plot(pbranchtab, labels=labels0,labeldirections=[horizontal,horizontal], color=colorlist[1..numelems(res)], transparency=0.5);
    logtrajfig:=Show2DLogTraj(file2,eps,nbcolumn,listnotread,name1,name2);
    figure:=plots[display](fig,logtrajfig);
    return figure;
end:



WriteVRepresentation:=proc(fileInput,fileOutput)
	#Write a file describing the V-representation of a given file describing polyhedra.
    local branches;
    with(PolyhedralSets);
    branches:=MakePolyhedra(fileInput);
    WriteVRepr(branches,fileOutput);
end:



WriteVRepr:=proc(branches,fileOutput)
	#Write a file describing the V-representation of a given file describing polyhedra, subprocess.
    local v, r, b, s1, s2, e, ve, re, i, j, k, sr;
    FileTools[Text][Open](fileOutput,append=true,overwrite=true);
    for i from 1 to numelems(branches) do
        b:=branches[i];
        sr:="\nRays:\n";
        v,r:=VerticesAndRays(b);
        s1:=cat("\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ file ",convert(i-1,string),"\n\nVertices:\n");
        FileTools[Text][WriteString](fileOutput,s1);
        for j from 1 to numelems(v) do
            s2:=cat(convert(v[j],string),"\n");
            FileTools[Text][WriteString](fileOutput,s2);
        od;
        FileTools[Text][WriteString](fileOutput,"\nEdges:\n");
        e:=Edges(b);
        for k from 1 to numelems(e) do
            ve,re:=VerticesAndRays(e[k]);
            if re=[] then
                s2:=cat(convert(ve,string),"\n");
                FileTools[Text][WriteString](fileOutput,s2);
            else
                sr:=cat(sr,substring(convert(ve,string),1..-2,string),", ",substring(convert(re,string),2..-1,string),"\n");
            end if;
        od;
        FileTools[Text][WriteString](fileOutput,sr);
    od;
    FileTools[Text][Close](fileOutput);
end:



TakeDimInter:=proc(ps1,ps2)
	#From two list ps1 and ps2 of polyhedra, return the dimension of the biggest intersection between one polyhedra of ps1 and one of ps2.
    local i, j, maxl, maxg, p;
    maxg:=-1;
    for i from 1 to numelems(ps1) do
        for j from 1 to numelems(ps2) do
            p:= ps1[i] intersect ps2[j];
            maxl:=Dimension(p);
            if maxl>maxg then
                maxg:=maxl;
            end if;
        od;
    od;
    return maxg;
end:



MakeListOfInter:=proc(ps1,ps2)
	#From two list ps1 and ps2 of polyhedra, return the list of intersections between ps1's polyhedra and ps2's polyhedra, in the form (index for ps1, index for ps2, polyhedron p1 in ps1, dimension of p1, polyhedron p2 in ps2, dimension of p2, instersected polyhedron p, dimension of p).
    local d, d1, d2, i, j, p, l;
    l:=[];
    for i from 1 to numelems(ps1) do
        for j from 1 to numelems(ps2) do
            p:= ps1[i] intersect ps2[j];
            d:= Dimension(p);
            d1:= Dimension(ps1[i]);
            d2:= Dimension(ps2[j]);
            if d>-1 then
                l:=[op(l),[i,j,ps1[i],d1,ps2[j],d2,p,d]];
            end if;
        od;
    od;
    return l;
end:

TakeMaxInter:=proc(ps1,ps2)
	#From two list ps1 and ps2 of polyhedra, return the list of intersections of maximal dimension between ps1's polyhedra and ps2's polyhedra, in the form (index for ps1, index for ps2, polyhedron p1 in ps1, dimension of p1, polyhedron p2 in ps2, dimension of p2, instersected polyhedron p, dimension of p).
    local d, d1, d2, i, j, p, l, dmax;
    l:=[];
    dmax:=-1;
    for i from 1 to numelems(ps1) do
        for j from 1 to numelems(ps2) do
            p:= ps1[i] intersect ps2[j];
            d:= Dimension(p);
            d1:= Dimension(ps1[i]);
            d2:= Dimension(ps2[j]);
            if d>-1 then
                if d=dmax then
                    l:=[op(l),[i,j,ps1[i],d1,ps2[j],d2,p,d]];
                else 
                    if d>dmax then
                        dmax:=d;
                        l:=[];
                        l:=[op(l),[i,j,ps1[i],d1,ps2[j],d2,p,d]];
                    end if;
                end if;
            end if;
        od;
    od;
    return l;
end:



MakeInterPolyhedra:=proc(interlist)
	#From a result given by procedures TakeMaxInter or MakeListOfInter, gives only instersected polyhedra.
	local c, l;
	l:=[];
	for c from 1 to numelems(interlist) do
		l:=[op(l),interlist[c][7]];
	od;
	return l;
end:



ShowMultipleTrace:=proc(eps,nbrcolumn,listnotread)
    #This is a work function, not a general one. Show trajectories of the columns to be read in the given datafile in log[eps] basis. nbrcolumn is the number of columns in the given datafile, and listnotread is a list of columns to be deleted from the data. It should remain 2 columns. 
    local i, d, dd, ddd, logtraj, A, C, filelist, stringfile;
    with(plots);
    filelist:=[];
    A:=[];
    for i from 1 to 9 do
    	stringfile:=cat("./00",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 10 to 99 do
    	stringfile:=cat("./0",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 100 to 375 do
    	stringfile:=cat("./",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 1 to 375 do
	    d:=readdata(filelist[i],nbrcolumn);
	    ddd:=LinearAlgebra[DeleteColumn](d,listnotread);
	    dd:=LinearAlgebra[DeleteRow](ddd,[1,2]);
	    logtraj:=map(log[eps],dd);
	    #B:=plots[textplot]([[logtraj[1,1],logtraj[1,2],"alpha",color="Red"],[logtraj[-1,1],logtraj[-1,2],"omega",color="Blue"]]);
	    A:=[op(A),plots[pointplot](logtraj,color="Black",symbolsize=2)];
	od;
    C:=plots[display](A,labels=["a3","a4"],labeldirections=[horizontal,horizontal]);
    return C;
end:



LogRoundTraj2D:=proc(file,eps,prec,nbrcolumn,listnotread,name1,name2)
    #Show trajectories of the columns to be read in the given datafile in a rounded log[eps] basis (if prec=1, it will give integers). nbrcolumn is the number of columns in the given datafile, and listnotread is a list of columns to be deleted from the data. It should remain 2 columns.
    local d, dd, ddd, logtraj, B, A, C, labels0;
    with(plots);
    d:=readdata(file,nbrcolumn);
    ddd:=LinearAlgebra[DeleteColumn](d,listnotread);
    dd:=LinearAlgebra[DeleteRow](ddd,[1]);
    logtraj:=map(x->(round(prec*log[eps](x))/prec),dd);
    labels0:=[name1,name2];
    #B:=plots[textplot]([[logtraj[1,1],logtraj[1,2],"alpha",color="Red"],[logtraj[-1,1],logtraj[-1,2],"omega",color="Blue"]]);
    #A:=plots[pointplot](logtraj,color="Black",symbolsize=10);
    A:=plots[pointplot](logtraj,color="Red",symbol=solidcircle,symbolsize=10,labels=labels0);
    C:=plots[display](A);#,B);
    return C;
end:



TraceAnalyser:=proc(file,eps,prec,params)
	#This is a work function, not a general one. From the given trace, give the list of tropical equilibrations "encountered" and the the number of points spent in it.
	local n, i, d, dd, ddd, logtrace, logtraj, logparams, vx1, vx2, vx3, vx4, vx5, pol1, pol2, pol3, pol4, pol5, pol6, pol7, orders, eq1, eq2, listeq, cpt;
	ddd:=readdata(file,6);
	dd:=LinearAlgebra[DeleteColumn](ddd,1);
	d:=LinearAlgebra[DeleteRow](dd,1);
	logtrace:=map(x->round(prec*log[eps](x))/prec,d);
	#print(logtrace);
	logparams:=map(x->round(prec*log[eps](x))/prec,params);
	#print(logparams);
	vx1:=LinearAlgebra[DeleteColumn](logtrace,[2,3,4,5]);
	vx2:=LinearAlgebra[DeleteColumn](logtrace,[1,3,4,5]);
	vx3:=LinearAlgebra[DeleteColumn](logtrace,[1,2,4,5]);
	vx4:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,5]);
	vx5:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,4]);
	pol1:=logparams[1]+~vx3;
	#k1x3
	#print(pol1);
	pol2:=logparams[2]+~vx1;
	#k2x1
	pol3:=logparams[3]+~vx2;
	#k3x2
	pol4:=logparams[4]+~vx2+~vx5;
	#k4x2x5
	pol5:=map(x->logparams[5],vx1);
	#k6
	#print(pol5);
	pol6:=logparams[6]+~vx3+~vx3+~vx4;
	#k9x3^2x4
	pol7:=logparams[7]+~vx4;
	#k10x4
	n:=LinearAlgebra[RowDimension](d);
	#orders:=[];
	listeq:=[];
	eq1:="";
	cpt:=0;
	for i from 1 to n do
		orders:=[min(pol1[i][1],pol2[i][1],pol3[i][1])-vx1[i][1], min(pol2[i][1],pol3[i][1],pol4[i][1])-vx2[i][1], min(pol1[i][1],pol6[i][1],pol7[i][1])-vx3[i][1], min(pol4[i][1],pol6[i][1],pol7[i][1])-vx4[i][1], min(pol5[i][1], pol4[i][1])-vx5[i][1]];
		#print(orders);
		eq2:=CheckEquilibration([pol1[i][1],pol2[i][1],pol3[i][1],pol4[i][1],pol5[i][1],pol6[i][1],pol7[i][1]],orders);
		#print(eq2);
		if (eq1<>eq2) then
			listeq:=[op(listeq),[eq1,cpt]];
			eq1:=eq2;
			cpt:=1;
		else
			cpt:=cpt+1;
		end if;
	od;
	listeq:=[op(listeq),[eq1,cpt]];
	return listeq;
end:



CheckEquilibration:=proc(pols,ords)
	#This is a work function, not a general one. From a point, check in which tropical equilibration we are.
	local b1, b2, b3, b4, b5, eqbtot, eqbx3, eqbx4, eqbx5, eqbx3x5, eqbx3x4, eqbx3x4x5, eqbx1x3x4x5, eqbx2x3x4x5;
	if (min(pols[1],pols[3]) = pols[2]) then
		b1:=true;
	else
		b1:=false;
	end if;
	if (min(pols[4],pols[3]) = pols[2]) then
		b2:=true;
	else
		b2:=false;
	end if;
	if (min(pols[6],pols[7]) = pols[1]) then
		b3:=true;
	else
		b3:=false;
	end if;
	if (min(pols[6],pols[7]) = pols[4]) then
		b4:=true;
	else
		b4:=false;
	end if;
	if (pols[4] = pols[5]) then
		b5:=true;
	else
		b5:=false;
	end if;
	if (b1 and b2 and b3 and b4 and b5) then
		eqbtot:=true;
	else
		eqbtot:=false;
	end if;
	if (b1 and b2 and b4 and b5 and (ords[3]>=max(ords[1],ords[2],ords[4],ords[5]))) then
		eqbx3:=true;
	else
		eqbx3:=false;
	end if;
	if (b1 and b2 and b3 and b5 and (ords[4]>=max(ords[1],ords[2],ords[3],ords[5]))) then
		eqbx4:=true;
	else
		eqbx4:=false;
	end if;
	if (b1 and b2 and b4 and b3 and (ords[5]>=max(ords[1],ords[2],ords[4],ords[3]))) then
		eqbx5:=true;
	else
		eqbx5:=false;
	end if;
	if (b1 and b2 and b5 and (min(ords[3],ords[4])>=max(ords[1],ords[2],ords[5]))) then
		eqbx3x4:=true;
	else
		eqbx3x4:=false;
	end if;
	if (b1 and b2 and b4 and (min(ords[3],ords[5])>=max(ords[1],ords[2],ords[4]))) then
		eqbx3x5:=true;
	else
		eqbx3x5:=false;
	end if;
	if (b1 and b2 and (min(ords[3],ords[4],ords[5])>=max(ords[1],ords[2]))) then
		eqbx3x4x5:=true;
	else
		eqbx3x4x5:=false;
	end if;
	if (b1 and (min(ords[2],ords[3],ords[4],ords[5])>=ords[1])) then
		eqbx2x3x4x5:=true;
	else
		eqbx2x3x4x5:=false;
	end if;
	if (b2 and (min(ords[1],ords[3],ords[4],ords[5])>=ords[2])) then
		eqbx1x3x4x5:=true;
	else
		eqbx1x3x4x5:=false;
	end if;
	if eqbtot then
		return "btot";
	else
		if eqbx3 then
			return "bx3";
		else
			if eqbx4 then
				return "bx4";
			else
				if eqbx5 then
					return "bx5";
				else
					if eqbx3x4 then
						return "bx3x4";
					else
						if eqbx3x5 then
							return "bx3x5";
						else
							if eqbx3x4x5 then
								return "bx3x4x5";
							else
								if eqbx2x3x4x5 then
									return "bx2x3x4x5";
								else
									if eqbx1x3x4x5 then
										return "bx1x3x4x5";
									else
										return "outside";
									end if;
								end if;
							end if;
						end if;
					end if;
				end if;
			end if;
		end if;
	end if;
end:



TraceAnalyser2:=proc(file,eps,prec,params)
	#This is a work function, not a general one. Obsolete, do not use it as the rounding is made later. From the given trace, give the list of tropical equilibrations "encountered" and the the number of points spent in it.
	local n, i, d, dd, ddd, logtrace, logtraj, logparams, vx1, vx2, vx3, vx4, vx5, pol1, pol2, pol3, pol4, pol5, pol6, pol7, orders, eq1, eq2, listeq, cpt;
	ddd:=readdata(file,6);
	dd:=LinearAlgebra[DeleteColumn](ddd,1);
	d:=LinearAlgebra[DeleteRow](dd,1);
	logtrace:=map(x->log[eps](x),d);
	#print(logtrace);
	logparams:=map(x->convert(log[eps](x),float),params);
	#print(logparams);
	vx1:=LinearAlgebra[DeleteColumn](logtrace,[2,3,4,5]);
	vx2:=LinearAlgebra[DeleteColumn](logtrace,[1,3,4,5]);
	vx3:=LinearAlgebra[DeleteColumn](logtrace,[1,2,4,5]);
	vx4:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,5]);
	vx5:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,4]);
	pol1:=logparams[1]+~vx3;
	#k1x3
	#print(pol1);
	pol2:=logparams[2]+~vx1;
	#k2x1
	#print(pol2);
	pol3:=logparams[3]+~vx2;
	#k3x2
	pol4:=logparams[4]+~vx2+~vx5;
	#k4x2x5
	pol5:=map(x->logparams[5],vx1);
	#k6
	#print(pol5);
	pol6:=logparams[6]+~vx3+~vx3+~vx4;
	#k9x3^2x4
	pol7:=logparams[7]+~vx4;
	#k10x4
	n:=LinearAlgebra[RowDimension](d);
	#orders:=[];
	listeq:=[];
	eq1:="";
	cpt:=0;
	for i from 1 to n do
		orders:=[min(pol1[i][1],pol2[i][1],pol3[i][1])-vx1[i][1], min(pol2[i][1],pol3[i][1],pol4[i][1])-vx2[i][1], min(pol1[i][1],pol6[i][1],pol7[i][1])-vx3[i][1], min(pol4[i][1],pol6[i][1],pol7[i][1])-vx4[i][1], min(pol5[i][1], pol4[i][1])-vx5[i][1]];
		eq2:=CheckEquilibration2([round(prec*pol1[i][1])/prec,round(prec*pol2[i][1])/prec,round(prec*pol3[i][1])/prec,round(prec*pol4[i][1])/prec,round(prec*pol5[i][1])/prec,round(prec*pol6[i][1])/prec,round(prec*pol7[i][1])/prec],orders);
		#print(eq2);
		if (eq1<>eq2) then
			listeq:=[op(listeq),[eq1,cpt]];
			eq1:=eq2;
			cpt:=1;
		else
			cpt:=cpt+1;
		end if;
	od;
	listeq:=[op(listeq),[eq1,cpt]];
	return listeq;
end:



CheckEquilibration2:=proc(pols,ords)
	#This is a work function, not a general one. Obsolete, do not use it as the rounding is made later. From a point, check in which tropical equilibration we are.
	local b1, b2, b3, b4, b5, eqbtot, eqbx3, eqbx4, eqbx5, eqbx3x5, eqbx3x4, eqbx3x4x5, eqbx1x3x4x5, eqbx2x3x4x5;
	if (min(pols[1],pols[3]) = pols[2]) then
		b1:=true;
	else
		b1:=false;
	end if;
	if (min(pols[4],pols[3]) = pols[2]) then
		b2:=true;
	else
		b2:=false;
	end if;
	if (min(pols[6],pols[7]) = pols[1]) then
		b3:=true;
	else
		b3:=false;
	end if;
	if (min(pols[6],pols[7]) = pols[4]) then
		b4:=true;
	else
		b4:=false;
	end if;
	if (pols[4] = pols[5]) then
		b5:=true;
	else
		b5:=false;
	end if;
	if (b1 and b2 and b3 and b4 and b5) then
		eqbtot:=true;
	else
		eqbtot:=false;
	end if;
	if (b1 and b2 and b4 and b5 and (ords[3]>=max(ords[1],ords[2],ords[4],ords[5]))) then
		eqbx3:=true;
	else
		eqbx3:=false;
	end if;
	if (b1 and b2 and b3 and b5 and (ords[4]>=max(ords[1],ords[2],ords[3],ords[5]))) then
		eqbx4:=true;
	else
		eqbx4:=false;
	end if;
	if (b1 and b2 and b4 and b3 and (ords[5]>=max(ords[1],ords[2],ords[4],ords[3]))) then
		eqbx5:=true;
	else
		eqbx5:=false;
	end if;
	if (b1 and b2 and b5 and (min(ords[3],ords[4])>=max(ords[1],ords[2],ords[5]))) then
		eqbx3x4:=true;
	else
		eqbx3x4:=false;
	end if;
	if (b1 and b2 and b4 and (min(ords[3],ords[5])>=max(ords[1],ords[2],ords[4]))) then
		eqbx3x5:=true;
	else
		eqbx3x5:=false;
	end if;
	if (b1 and b2 and (min(ords[3],ords[4],ords[5])>=max(ords[1],ords[2]))) then
		eqbx3x4x5:=true;
	else
		eqbx3x4x5:=false;
	end if;
	if (b1 and (min(ords[2],ords[3],ords[4],ords[5])>=ords[1])) then
		eqbx2x3x4x5:=true;
	else
		eqbx2x3x4x5:=false;
	end if;
	if (b2 and (min(ords[1],ords[3],ords[4],ords[5])>=ords[2])) then
		eqbx1x3x4x5:=true;
	else
		eqbx1x3x4x5:=false;
	end if;
	if eqbtot then
		return "btot";
	else
		if eqbx3 then
			return "bx3";
		else
			if eqbx4 then
				return "bx4";
			else
				if eqbx5 then
					return "bx5";
				else
					if eqbx3x4 then
						return "bx3x4";
					else
						if eqbx3x5 then
							return "bx3x5";
						else
							if eqbx3x4x5 then
								return "bx3x4x5";
							else
								if eqbx2x3x4x5 then
									return "bx2x3x4x5";
								else
									if eqbx1x3x4x5 then
										return "bx1x3x4x5";
									else
										return "outside";
									end if;
								end if;
							end if;
						end if;
					end if;
				end if;
			end if;
		end if;
	end if;
end:



ShowMultipleTrace2:=proc(eps,nbrcolumn,listnotread)
    #This is a work function, not a general one. Show trajectories of the columns to be read in the given datafile in log[eps] basis. nbrcolumn is the number of columns in the given datafile, and listnotread is a list of columns to be deleted from the data. It should remain 3 columns. 
    local j, n, l, i, d, dd, ddd, logtraj, A, C, filelist, stringfile;
    with(plots);
    filelist:=[];
    A:=[];
    for i from 1 to 9 do
    	stringfile:=cat("./00",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 10 to 99 do
    	stringfile:=cat("./0",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 100 to 375 do
    	stringfile:=cat("./",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 2 to 375 do
	    d:=readdata(filelist[i],nbrcolumn);
		ddd:=LinearAlgebra[DeleteColumn](d,listnotread);
    	#n:=LinearAlgebra[RowDimension](ddd);
    	#j:=1;
    	#l:=[];
    	#while j+9<n do
	    #	l:=[op(l),j,j+1,j+3,j+4,j+5,j+6,j+7,j+8,j+9];
	    #	j:=j+10;
	    #od;
    	dd:=LinearAlgebra[DeleteRow](ddd,1);
    	logtraj:=map(log[eps],dd);
	    #B:=plots[textplot]([[logtraj[1,1],logtraj[1,2],"alpha",color="Red"],[logtraj[-1,1],logtraj[-1,2],"omega",color="Blue"]]);
	    A:=[op(A),plots[pointplot3d](logtraj,color="Black",symbol=solidsphere,symbolsize=2,size=[0.8,0.8],orientation=[109,13,-47])];
	od;
	d:=readdata(filelist[1],nbrcolumn);
	ddd:=LinearAlgebra[DeleteColumn](d,listnotread);
    dd:=LinearAlgebra[DeleteRow](ddd,1);
    logtraj:=map(log[eps],dd);
	A:=[op(A),plots[pointplot3d](logtraj,color="Red",symbolsize=3,size=[0.8,0.8],orientation=[109,13,-47])];
    C:=plots[display](A,labels=["a3","a4","a5"],labeldirections=[horizontal,horizontal,horizontal],size=[0.8,0.8],orientation=[109,13,-47]);
    return C;
end:



Show2Dtrace:=proc(file,eps,nbrcolumn,listnotread,mcolor:="Black")
    # Show trajectories of the columns to be read in the given datafile in log[eps] basis. nbrcolumn is the number of columns in the given datafile, and listnotread is a list of columns to be deleted from the data. It should remain 2 columns.
    local d, dd, ddd, logtraj, B, A, C;
    with(plots);
    d:=readdata(file,nbrcolumn);
    ddd:=LinearAlgebra[DeleteColumn](d,listnotread);
    dd:=LinearAlgebra[DeleteRow](ddd,1);
    logtraj:=map(log[eps],dd);
    #B:=plots[textplot]([[logtraj[1,1],logtraj[1,2],"alpha",color="Red"],[logtraj[-1,1],logtraj[-1,2],"omega",color="Blue"]]);
    A:=plots[pointplot](logtraj,color=mcolor,symbolsize=5);
    #C:=plots[display](A,B);
    return A;
end:


CheckEquilibrationAndBranchesd1e11:=proc(pols,ords,species)
	#This is a work function, not a general one. A more detailed analysis of a point: gives the tropical equilibration where the tropicalized point lives and the polytope in the format nnnnnebb when n is the species (0 for x_i fast, i for x_i slow) e is an error code (0 for fine, 1 for error, except in the case 12345100 which means 0 values in the trace), bb is the FILE bb in the given solution by smtcutpartialpol.py.
	local b1, b2, b3, b4, b5, eqbtot, eqbx3, eqbx4, eqbx5, eqbx3x5, eqbx3x4, eqbx3x4x5, eqbx1x3x4x5, eqbx2x3x4x5, p1p2, p2p3, p2p4, p1p6, p1p7, p4p6, p4p7, s1, s2, s3, s4, x1p1, x1p2, x1p3, x2p2, x2p3, x2p4, x3p1, x3p6, x3p7, x4p4, x4p6, x4p7, x5p4, x5p5, p4p5, iszero;
	p1p2:=false;
	p2p3:=false;
	p2p4:=false;
	p1p6:=false;
	p1p7:=false;
	p4p6:=false;
	p4p7:=false;
	p4p5:=false;
	x1p1:=false;
	x1p2:=false;
	x1p3:=false;
	x2p2:=false;
	x2p3:=false;
	x2p4:=false;
	x3p1:=false;
	x3p6:=false;
	x3p7:=false;
	x4p4:=false;
	x4p6:=false;
	x4p7:=false;
	x5p4:=false;
	x5p5:=false;
	s1:=false;
	s2:=false;
	s3:=false;
	s4:=false;
	#print(species);
	if (evalb(species[1]=infinity) or evalb(species[2]=infinity) or evalb(species[3]=infinity) or evalb(species[4]=infinity) or evalb(species[5]=infinity)) then
		#print(y);
		return "12345100";
	end if;
	#if (evalb(ords[1]=infinity) or evalb(ords[2]=infinity) or evalb(ords[3]=infinity) or evalb(ords[4]=infinity) or evalb(ords[5]=infinity)) then
	#	return "12345100";
	#end if;
	if (min(pols[1],pols[3]) = pols[2]) then
		b1:=true;
		if (min(pols[1],pols[3]) = pols[1]) then
			p1p2:=true;
		else
			p2p3:=true;
		end if;
	else
		b1:=false;
	end if;
	if (min(pols[4],pols[3]) = pols[2]) then
		b2:=true;
		if (min(pols[4],pols[3]) = pols[4]) then
			p2p4:=true;
		else
			p2p3:=true;
		end if;
	else
		b2:=false;
	end if;
	if (min(pols[6],pols[7]) = pols[1]) then
		b3:=true;
		if (min(pols[6],pols[7]) = pols[6]) then
			p1p6:=true;
		else
			p1p7:=true;
		end if;
	else
		b3:=false;
	end if;
	if (min(pols[6],pols[7]) = pols[4]) then
		b4:=true;
		if (min(pols[6],pols[7]) = pols[6]) then
			p4p6:=true;
		else
			p4p7:=true;
		end if;
	else
		b4:=false;
	end if;
	if (pols[4] = pols[5]) then
		b5:=true;
		p4p5:=true;
	else
		b5:=false;
	end if;
	if (pols[1]<=pols[2] and pols[1]<=pols[3]) then
		x1p1:=true;
	end if;
	if (pols[2]<=pols[1] and pols[2]<=pols[3]) then
		x1p2:=true;
	end if;
	if (pols[3]<=pols[2] and pols[3]<=pols[1]) then
		x1p3:=true;
	end if;
	if (pols[2]<=pols[3] and pols[2]<=pols[4]) then
		x2p2:=true;
	end if;
	if (pols[3]<=pols[2] and pols[3]<=pols[4]) then
		x2p3:=true;
	end if;
	if (pols[4]<=pols[2] and pols[4]<=pols[3]) then
		x2p4:=true;
	end if;
	if (pols[1]<=pols[7] and pols[1]<=pols[6]) then
		x3p1:=true;
	end if;
	if (pols[6]<=pols[7] and pols[6]<=pols[1]) then
		x3p6:=true;
	end if;
	if (pols[7]<=pols[1] and pols[7]<=pols[6]) then
		x3p7:=true;
	end if;
	if (pols[4]<=pols[6] and pols[4]<=pols[7]) then
		x4p4:=true;
	end if;
	if (pols[6]<=pols[4] and pols[6]<=pols[7]) then
		x4p6:=true;
	end if;
	if (pols[7]<=pols[4] and pols[7]<=pols[6]) then
		x4p7:=true;
	end if;
	if (pols[4]<=pols[5]) then
		x5p4:=true;
	end if;
	if (pols[5]<=pols[4]) then
		x5p5:=true;
	end if;
	if (species[1]<=species[2] and species[1]<=species[3] and species[1]<=species[4]) then
		s1:=true;
	end if;
	if (species[2]<=species[1] and species[2]<=species[3] and species[2]<=species[4]) then
		s2:=true;
	end if;
	if (species[3]<=species[1] and species[3]<=species[2] and species[3]<=species[4]) then
		s3:=true;
	end if;
	if (species[4]<=species[1] and species[4]<=species[2] and species[4]<=species[3]) then
		s4:=true;
	end if;
	if (b1 and b2 and b3 and b4 and b5) then
		eqbtot:=true;
	else
		eqbtot:=false;
	end if;
	if (b1 and b2 and b4 and b5 and (ords[3]>=max(ords[1],ords[2],ords[4],ords[5]))) then
		eqbx3:=true;
	else
		eqbx3:=false;
	end if;
	if (b1 and b2 and b3 and b5 and (ords[4]>=max(ords[1],ords[2],ords[3],ords[5]))) then
		eqbx4:=true;
	else
		eqbx4:=false;
	end if;
	if (b1 and b2 and b4 and b3 and (ords[5]>=max(ords[1],ords[2],ords[4],ords[3]))) then
		eqbx5:=true;
	else
		eqbx5:=false;
	end if;
	if (b1 and b2 and b5 and (min(ords[3],ords[4])>=max(ords[1],ords[2],ords[5]))) then
		eqbx3x4:=true;
	else
		eqbx3x4:=false;
	end if;
	if (b1 and b2 and b4 and (min(ords[3],ords[5])>=max(ords[1],ords[2],ords[4]))) then
		eqbx3x5:=true;
	else
		eqbx3x5:=false;
	end if;
	if (b1 and b2 and (min(ords[3],ords[4],ords[5])>=max(ords[1],ords[2]))) then
		eqbx3x4x5:=true;
	else
		eqbx3x4x5:=false;
	end if;
	if (b1 and (min(ords[2],ords[3],ords[4],ords[5])>=ords[1])) then
		eqbx2x3x4x5:=true;
	else
		eqbx2x3x4x5:=false;
	end if;
	if (b2 and (min(ords[1],ords[3],ords[4],ords[5])>=ords[2])) then
		eqbx1x3x4x5:=true;
	else
		eqbx1x3x4x5:=false;
	end if;

	if eqbtot then
		if (p2p3 and p2p3 and p1p6 and p4p6 and p4p5 and s4) then
			return "00000000";# "btot_00";
		elif (p1p2 and p2p4 and p1p6 and p4p6 and p4p5 and s4) then
			return "00000001";# "btot_01";
		else
			return "00000100";# "btot_??";
		end if;
	elif eqbx3 then
		if (p2p3 and p2p3 and x3p1 and p4p6 and p4p5 and s2) then
			return "00300000";# "bx3_00";
		elif (p2p3 and p2p3 and x3p1 and p4p6 and p4p5 and s3) then
			return "00300001";# "bx3_01";
		else
			return "00300100";# "bx3_??";
		end if;
	elif eqbx4 then
		if (p2p3 and p2p3 and p1p7 and x4p4 and p4p5 and s2) then
			return "00040000";# "bx4_00";
		elif (p2p3 and p2p3 and p1p7 and x4p6 and p4p5 and s4) then
			return "00040001";# "bx4_01";
		elif (p2p3 and p2p3 and p1p6 and x4p6 and p4p5 and s2) then
			return "00040002";# "bx4_02";
		else
			return "00040100";# "bx4_??";
		end if;
	elif eqbx5 then
		if (p2p3 and p2p3 and p1p6 and p4p7 and x5p4 and s4) then
			return "00005000";# "bx5_00";
		elif (p1p2 and p2p4 and p1p6 and p4p7 and x5p4 and s4) then
			return "00005001";# "bx5_01";
		elif (p1p2 and p2p4 and p1p6 and p4p6 and x5p4 and s3) then
			return "00005002";# "bx5_02";
		elif (p2p3 and p2p3 and p1p6 and p4p6 and x5p4 and s3) then
			return "00005003";# "bx5_03";
		else
			return "00005100";# "bx5_??";
		end if;
	elif eqbx3x4 then
		if (p2p3 and p2p3 and x3p6 and x4p6 and p4p5 and s2) then
			return "00340000";# "bx3x4_00";
		elif (p2p3 and p2p3 and x3p6 and x4p6 and p4p5 and s4) then
			return "00340001";# "bx3x4_01";
		elif (p2p3 and p2p3 and x3p7 and x4p7 and p4p5 and s4) then
			return "00340002";# "bx3x4_02";
		elif (p2p3 and p2p3 and x3p7 and x4p4 and p4p5 and s2) then
			return "00340003";# "bx3x4_03";
		elif (p2p3 and p2p3 and x3p1 and x4p4 and p4p5 and s2) then
			return "00340004";# "bx3x4_04";
		elif (p2p3 and p2p3 and x3p1 and x4p6 and p4p5 and s2) then
			return "00340005";# "bx3x4_05";
		else
			return "00340100";# "bx3x4_??";
		end if;
	elif eqbx3x5 then
		if (p1p2 and p2p4 and x3p1 and p4p6 and x5p4 and s3) then
			return "00305000";# "bx3x5_00";
		elif (p2p3 and p2p3 and x3p1 and p4p6 and x5p4 and s3) then
			return "00305001";# "bx3x5_01";
		elif (p2p3 and p2p3 and x3p1 and p4p6 and x5p5 and s3) then
			return "00305002";# "bx3x5_02";
		elif (p2p3 and p2p3 and x3p6 and p4p6 and x5p4 and s3) then
			return "00305003";# "bx3x5_03";
		else
			return "00305100";# "bx3x5_??";
		end if;
	elif eqbx3x4x5 then
		if (p2p3 and p2p3 and x3p1 and x4p4 and x5p5 and s2) then
			return "00345000";# "bx3x4x5_00";
		elif (p2p3 and p2p3 and x3p1 and x4p4 and x5p4 and s2) then
			return "00345001";# "bx3x4x5_01";
		elif (p2p3 and p2p3 and x3p7 and x4p4 and x5p4 and s2) then
			return "00345002";# "bx3x4x5_02";
		elif (p2p3 and p2p3 and x3p7 and x4p4 and x5p5 and s2) then
			return "00345003";# "bx3x4x5_03";
		elif (p2p3 and p2p3 and x3p7 and x4p7 and x5p5 and s2) then
			return "00345004";# "bx3x4x5_04";
		elif (p2p3 and p2p3 and x3p1 and x4p6 and x5p5 and s2) then
			return "00345005";# "bx3x4x5_05";
		elif (p2p3 and p2p3 and x3p1 and x4p7 and x5p5 and s2) then
			return "00345006";# "bx3x4x5_06";
		elif (p2p3 and p2p3 and x3p6 and x4p4 and x5p4 and s2) then
			return "00345007";# "bx3x4x5_07";
		elif (p2p3 and p2p3 and x3p6 and x4p4 and x5p4 and s3) then
			return "00345008";# "bx3x4x5_08";
		elif (p2p3 and p2p3 and x3p6 and x4p4 and x5p4 and s4) then
			return "00345009";# "bx3x4x5_09";
		elif (p2p3 and p2p3 and x3p6 and x4p6 and x5p5 and s2) then
			return "00345010";# "bx3x4x5_10";
		elif (p2p3 and p2p3 and x3p6 and x4p6 and x5p4 and s2) then
			return "00345011";# "bx3x4x5_11";
		elif (p2p3 and p2p3 and x3p1 and x4p6 and x5p4 and s2) then
			return "00345012";# "bx3x4x5_12";
		elif (p1p2 and p2p4 and x3p1 and x4p4 and x5p4 and s3) then
			return "00345013";# "bx3x4x5_13";
		elif (p1p2 and p2p4 and x3p7 and x4p7 and x5p5 and s4) then
			return "00345014";# "bx3x4x5_14";
		elif (p1p2 and p2p4 and x3p6 and x4p6 and x5p4 and s3) then
			return "00345015";# "bx3x4x5_15";
		elif (p1p2 and p2p4 and x3p6 and x4p6 and x5p4 and s4) then
			return "00345016";# "bx3x4x5_16";
		elif (p2p3 and p2p3 and x3p6 and x4p6 and x5p4 and s4) then
			return "00345017";# "bx3x4x5_17";
		elif (p2p3 and p2p3 and x3p7 and x4p7 and x5p5 and s4) then
			return "00345018";# "bx3x4x5_18";
		elif (p2p3 and p2p3 and x3p6 and x3p6 and x5p5 and s4) then
			return "00345019";# "bx3x4x5_19";
		elif (p2p3 and p2p3 and x3p7 and x4p4 and x5p4 and s4) then
			return "00345020";# "bx3x4x5_20";
		elif (p2p3 and p2p3 and x3p6 and x4p6 and x5p5 and s3) then
			return "00345021";# "bx3x4x5_21";
		elif (p2p3 and p2p3 and x3p1 and x4p6 and x5p5 and s3) then
			return "00345022";# "bx3x4x5_22";
		elif (p2p3 and p2p3 and x3p1 and x4p4 and x5p5 and s3) then
			return "00345023";# "bx3x4x5_23";
		elif (p2p3 and p2p3 and x3p1 and x4p4 and x5p4 and s3) then
			return "00345024";# "bx3x4x5_24";
		elif (p2p3 and p2p3 and x3p1 and x4p6 and x5p4 and s3) then
			return "00345025";# "bx3x4x5_25";
		elif (p2p3 and p2p3 and x3p6 and x4p6 and x5p4 and s3) then
			return "00345026";# "bx3x4x5_26";
		else
			return "00345100";# "bx3x4x5_??";
		end if;
	elif eqbx1x3x4x5 then
		if (x1p2 and p2p4 and x3p7 and x4p4 and x5p4 and s1) then
			return "10345000";# "bx1x3x4x5_00";
		elif (x1p2 and p2p4 and x3p7 and x4p4 and x5p4 and s4) then
			return "10345001";# "bx1x3x4x5_01";
		elif (x1p2 and p2p4 and x3p1 and x4p4 and x5p4 and s3) then
			return "10345002";# "bx1x3x4x5_02";
		elif (x1p2 and p2p4 and x3p6 and x4p4 and x5p4 and s4) then
			return "10345003";# "bx1x3x4x5_03";
		elif (x1p2 and p2p4 and x3p6 and x4p6 and x5p4 and s4) then
			return "10345004";# "bx1x3x4x5_04";
		elif (x1p1 and p2p4 and x3p6 and x4p6 and x5p4 and s4) then
			return "10345005";# "bx1x3x4x5_05";
		elif (x1p1 and p2p4 and x3p6 and x4p6 and x5p5 and s4) then
			return "10345006";# "bx1x2x3x4_06";
		elif (x1p1 and p2p4 and x3p6 and x4p6 and x5p5 and s3) then
			return "10345007";# "bx1x3x4x5_07";
		elif (x1p1 and p2p4 and x3p6 and x4p6 and x5p5 and s3) then
			return "10345008";# "bx1x3x4x5_08";
		elif (x1p1 and p2p4 and x3p1 and x4p6 and x5p5 and s3) then
			return "10345009";# "bx1x3x4x5_09";
		elif (x1p1 and p2p4 and x3p1 and x4p4 and x5p4 and s3) then
			return "10345010";# "bx1x3x4x5_10";
		elif (x1p1 and p2p4 and x3p1 and x4p6 and x5p4 and s3) then
			return "10345011";# "bx1x3x4x5_11";
		elif (x1p1 and p2p4 and x3p6 and x4p6 and x5p4 and s3) then
			return "10345012";# "bx1x3x4x5_12";
		elif (x1p2 and p2p4 and x3p6 and x4p4 and x5p4 and s3) then
			return "10345013";# "bx1x3x4x5_13";
		elif (x1p2 and p2p4 and x3p6 and x4p6 and x5p4 and s3) then
			return "10345014";# "bx1x3x4x5_14";
		elif (x1p2 and p2p4 and x3p6 and x4p4 and x5p4 and s1) then
			return "10345015";# "bx1x3x4x5_15";
		elif (x1p2 and p2p4 and x3p1 and x4p4 and x5p4 and s1) then
			return "10345016";# "bx1x3x4x5_16";
		elif (x1p2 and p2p4 and x3p7 and x4p7 and x5p5 and s4) then
			return "10345017";# "bx1x3x4x5_17";
		elif (x1p1 and p2p4 and x3p7 and x4p7 and x5p5 and s4) then
			return "10345018";# "bx1x3x4x5_18";
		else
			return "10345100";# "bx1x3x4x5_??";
		end if;
	elif eqbx2x3x4x5 then
		if (p1p2 and x2p2 and x3p6 and x4p6 and x5p5 and s4) then
			return "02345000";# "bx2x3x4x5_00";
		elif (p1p2 and x2p2 and x3p6 and x4p6 and x5p4 and s4) then
			return "02345001";# "bx2x3x4x5_01";
		elif (p1p2 and x2p4 and x3p6 and x4p6 and x5p4 and s4) then
			return "02345002";# "bx2x3x4x5_02";
		elif (p1p2 and x2p4 and x3p6 and x4p4 and x5p4 and s4) then
			return "02345003";# "bx2x3x4x5_03";
		elif (p1p2 and x2p4 and x3p6 and x4p6 and x5p4 and s3) then
			return "02345004";# "bx2x3x4x5_04";
		elif (p1p2 and x2p4 and x3p1 and x4p4 and x5p4 and s3) then
			return "02345005";# "bx2x3x4x5_05";
		elif (p1p2 and x2p4 and x3p6 and x4p4 and x5p4 and s3) then
			return "02345006";# "bx2x3x4x5_06";
		elif (p1p2 and x2p2 and x3p6 and x4p6 and x5p5 and s3) then
			return "02345007";# "bx2x3x4x5_07";
		elif (p1p2 and x2p2 and x3p1 and x4p6 and x5p5 and s3) then
			return "02345008";# "bx2x3x4x5_08";
		elif (p1p2 and x2p2 and x3p1 and x4p4 and x5p5 and s3) then
			return "02345009";# "bx2x3x4x5_09";
		elif (p1p2 and x2p2 and x3p1 and x4p4 and x5p4 and s3) then
			return "02345010";# "bx2x3x4x5_10";
		elif (p1p2 and x2p2 and x3p6 and x4p6 and x5p4 and s3) then
			return "02345011";# "bx2x3x4x5_11";
		elif (p1p2 and x2p2 and x3p1 and x4p6 and x5p4 and s3) then
			return "02345012";# "bx2x3x4x5_12";
		elif (p1p2 and x2p4 and x3p7 and x4p4 and x5p4 and s4) then
			return "02345013";# "bx2x3x4x5_13";
		elif (p1p2 and x2p4 and x3p7 and x4p7 and x5p5 and s4) then
			return "02345014";# "bx2x3x4x5_14";
		elif (p1p2 and x2p2 and x3p7 and x4p7 and x5p5 and s4) then
			return "02345015";# "bx2x3x4x5_15";
		elif (p2p3 and x2p2 and x3p6 and x4p6 and x5p5 and s4) then
			return "02345016";# "bx2x3x4x5_16";
		elif (p2p3 and x2p2 and x3p6 and x4p6 and x5p4 and s4) then
			return "02345017";# "bx2x3x4x5_17";
		elif (p2p3 and x2p2 and x3p6 and x4p4 and x5p4 and s4) then
			return "02345018";# "bx2x3x4x5_18";
		elif (p2p3 and x2p2 and x3p6 and x4p4 and x5p4 and s2) then
			return "02345019";# "bx2x3x4x5_19";
		elif (p2p3 and x2p2 and x3p7 and x4p4 and x5p4 and s4) then
			return "02345020";# "bx2x3x4x5_20";
		elif (p2p3 and x2p4 and x3p7 and x4p4 and x5p4 and s4) then
			return "02345021";# "bx2x3x4x5_21";
		elif (p2p3 and x2p4 and x3p6 and x4p4 and x5p4 and s4) then
			return "02345022";# "bx2x3x4x5_22";
		elif (p2p3 and x2p4 and x3p6 and x4p6 and x5p4 and s4) then
			return "02345023";# "bx2x3x4x5_23";
		elif (p2p3 and x2p2 and x3p7 and x4p7 and x5p5 and s4) then
			return "02345024";# "bx2x3x4x5_24";
		elif (p2p3 and x2p4 and x3p7 and x4p7 and x5p5 and s4) then
			return "02345025";# "bx2x3x4x5_25";
		elif (p2p3 and x2p4 and x3p6 and x4p4 and x5p4 and s3) then
			return "02345026";# "bx2x3x4x5_26";
		elif (p2p3 and x2p4 and x3p1 and x4p4 and x5p4 and s3) then
			return "02345027";# "bx2x3x4x5_27";
		elif (p2p3 and x2p4 and x3p6 and x4p6 and x5p4 and s3) then
			return "02345028";# "bx2x3x4x5_28";
		elif (p2p3 and x2p2 and x3p6 and x4p6 and x5p4 and s3) then
			return "02345029";# "bx2x3x4x5_29";
		elif (p2p3 and x2p2 and x3p1 and x4p6 and x5p4 and s3) then
			return "02345030";# "bx2x3x4x5_30";
		elif (p2p3 and x2p2 and x3p6 and x4p6 and x5p5 and s3) then
			return "02345031";# "bx2x3x4x5_31";
		elif (p2p3 and x2p2 and x3p1 and x4p6 and x5p5 and s3) then
			return "02345032";# "bx2x3x4x5_32";
		elif (p2p3 and x2p2 and x3p1 and x4p4 and x5p4 and s3) then
			return "02345033";# "bx2x3x4x5_33";
		elif (p2p3 and x2p2 and x3p6 and x4p4 and x5p4 and s3) then
			return "02345034";# "bx2x3x4x5_34";
		elif (p2p3 and x2p2 and x3p1 and x4p4 and x5p5 and s2) then
			return "02345035";# "bx2x3x4x5_35";
		elif (p2p3 and x2p2 and x3p1 and x4p4 and x5p5 and s3) then
			return "02345036";# "bx2x3x4x5_36";
		elif (p2p3 and x2p4 and x3p6 and x4p4 and x5p4 and s2) then
			return "02345037";# "bx2x3x4x5_37";
		elif (p2p3 and x2p4 and x3p1 and x4p4 and x5p4 and s2) then
			return "02345038";# "bx2x3x4x5_38";
		elif (p2p3 and x2p4 and x3p7 and x4p4 and x5p4 and s2) then
			return "02345039";# "bx2x3x4x5_39";
		elif (p2p3 and x2p2 and x3p7 and x4p4 and x5p4 and s2) then
			return "02345040";# "bx2x3x4x5_40";
		elif (p2p3 and x2p2 and x3p7 and x4p4 and x5p5 and s2) then
			return "02345041";# "bx2x3x4x5_41";
		elif (p2p3 and x2p2 and x3p7 and x4p7 and x5p5 and s2) then
			return "02345042";# "bx2x3x4x5_42";
		elif (p2p3 and x2p2 and x3p1 and x4p7 and x5p5 and s2) then
			return "02345043";# "bx2x3x4x5_43";
		elif (p2p3 and x2p2 and x3p1 and x4p4 and x5p4 and s2) then
			return "02345044";# "bx2x3x4x5_44";
		elif (p2p3 and x2p2 and x3p1 and x4p6 and x5p4 and s2) then
			return "02345045";# "bx2x3x4x5_45";
		elif (p2p3 and x2p2 and x3p6 and x4p6 and x5p4 and s2) then
			return "02345046";# "bx2x3x4x5_46";
		elif (p2p3 and x2p2 and x3p6 and x4p6 and x5p5 and s2) then
			return "02345047";# "bx2x3x4x5_47";
		elif (p2p3 and x2p2 and x3p1 and x4p6 and x5p5 and s2) then
			return "02345048";# "bx2x3x4x5_48";
		else
			return "02345100";# "bx2x3x4x5_??";
		end if;
	else
		return "12345000";# "outside";
	end if;
end:

TraceAnalyserd1e11:=proc(file,eps,prec,params)
	#This is a work function, not a general one. A more detailed analysis of the trace, gives tropical equilibrations and the polytope in the format nnnnnebb when n is the species (0 for x_i fast, i for x_i slow) e is an error code (0 for fine, 1 for error, except in the case 12345100 which means 0 values in the trace), bb is the FILE bb in the given solution by smtcutpartialpol.py.
	local n, i, d, dd, ddd, logtrace, logtraj, logparams, vx1, vx2, vx3, vx4, vx5, pol1, pol2, pol3, pol4, pol5, pol6, pol7, orders, eq1, eq2, listeq, cpt, o1, o2, o3, o4, o5;
	ddd:=readdata(file,6);
	d:=LinearAlgebra[DeleteColumn](ddd,1);
	#d:=LinearAlgebra[DeleteRow](dd,1);
	logtrace:=map(x->round(prec*log[eps](abs(x)))/prec,d);
	#print(logtrace);
	logparams:=map(x->round(prec*log[eps](abs(x)))/prec,params);
	#print(logparams);
	vx1:=LinearAlgebra[DeleteColumn](logtrace,[2,3,4,5]);
	vx2:=LinearAlgebra[DeleteColumn](logtrace,[1,3,4,5]);
	vx3:=LinearAlgebra[DeleteColumn](logtrace,[1,2,4,5]);
	vx4:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,5]);
	vx5:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,4]);
	pol1:=logparams[1]+~vx3;
	#k1x3
	#print(pol1);
	pol2:=logparams[2]+~vx1;
	#k2x1
	pol3:=logparams[3]+~vx2;
	#k3x2
	pol4:=logparams[4]+~vx2+~vx5;
	#k4x2x5
	pol5:=map(x->logparams[5],vx1);
	#k6
	#print(pol5);
	pol6:=logparams[6]+~vx3+~vx3+~vx4;
	#k9x3^2x4
	pol7:=logparams[7]+~vx4;
	#k10x4
	n:=LinearAlgebra[RowDimension](d);
	#orders:=[];
	listeq:=[];
	eq1:="";
	cpt:=0;
	if (evalb(vx1[1][1]=infinity)) then 
		o1:=infinity;
	else
		o1:=min(pol1[1][1],pol2[1][1],pol3[1][1])-vx1[1][1];
	end if;
	if (evalb(vx2[1][1]=infinity)) then 
		o2:=infinity;
	else
		o2:=min(pol2[1][1],pol3[1][1],pol4[1][1])-vx2[1][1];
	end if;
	if (evalb(vx3[1][1]=infinity)) then 
		o3:=infinity;
	else
		o3:=min(pol1[1][1],pol6[1][1],pol7[1][1])-vx3[1][1];
	end if;
	if (evalb(vx4[1][1]=infinity)) then 
		o4:=infinity;
	else
		o4:=min(pol4[1][1],pol6[1][1],pol7[1][1])-vx4[1][1];
	end if;
	if (evalb(vx5[1][1]=infinity)) then 
		o5:=infinity;
	else
		o5:=min(pol5[1][1], pol4[1][1])-vx5[1][1];
	end if;
	orders:=[o1, o2, o3, o4, o5];
	#print(orders);
	eq1:=CheckEquilibrationAndBranchesd1e11([pol1[1][1],pol2[1][1],pol3[1][1],pol4[1][1],pol5[1][1],pol6[1][1],pol7[1][1]],orders,[vx1[1][1],vx2[1][1],vx3[1][1],vx4[1][1],vx5[1][1]]);
	cpt:=1;
	#print(eq2);
	for i from 2 to n do
		if (evalb(vx1[i][1]=infinity)) then 
			o1:=infinity;
		else
			o1:=min(pol1[i][1],pol2[i][1],pol3[i][1])-vx1[i][1];
		end if;
		if (evalb(vx2[i][1]=infinity)) then 
			o2:=infinity;
		else
			o2:=min(pol2[i][1],pol3[i][1],pol4[i][1])-vx2[i][1];
		end if;
		if (evalb(vx3[i][1]=infinity)) then 
			o3:=infinity;
		else
			o3:=min(pol1[i][1],pol6[i][1],pol7[i][1])-vx3[i][1];
		end if;
		if (evalb(vx4[i][1]=infinity)) then 
			o4:=infinity;
		else
			o4:=min(pol4[i][1],pol6[i][1],pol7[i][1])-vx4[i][1];
		end if;
		if (evalb(vx5[i][1]=infinity)) then 
			o5:=infinity;
		else
			o5:=min(pol5[i][1], pol4[i][1])-vx5[i][1];
		end if;
		orders:=[o1, o2, o3, o4, o5];
		#print(orders);
		eq2:=CheckEquilibrationAndBranchesd1e11([pol1[i][1],pol2[i][1],pol3[i][1],pol4[i][1],pol5[i][1],pol6[i][1],pol7[i][1]],orders,[vx1[i][1],vx2[i][1],vx3[i][1],vx4[i][1],vx5[i][1]]);
		#print(eq2);
		if (eq1<>eq2) then
			listeq:=[op(listeq),[eq1,cpt]];
			eq1:=eq2;
			cpt:=1;
		else
			cpt:=cpt+1;
		end if;
	od;
	listeq:=[op(listeq),[eq1,cpt]];
	return listeq;
end:

TraceAnalyser2d1e11:=proc(file,eps,prec,params)
	#This is a work function, not a general one. Obsolete, do not use it as the rounded is made later. A more detailed analysis of the trace, gives tropical equilibrations and the polytope in the format nnnnnebb when n is the species (0 for x_i fast, i for x_i slow) e is an error code (0 for fine, 1 for error, except in the case 12345100 which means 0 values in the trace), bb is the FILE bb in the given solution by smtcutpartialpol.py.
	local n, i, d, dd, ddd, logtrace, logtraj, logparams, vx1, vx2, vx3, vx4, vx5, pol1, pol2, pol3, pol4, pol5, pol6, pol7, orders, eq1, eq2, listeq, cpt;
	ddd:=readdata(file,6);
	dd:=LinearAlgebra[DeleteColumn](ddd,1);
	d:=LinearAlgebra[DeleteRow](dd,1);
	logtrace:=map(x->log[eps](abs(x)),d);
	#print(logtrace);
	logparams:=map(x->convert(log[eps](abs(x)),float),params);
	#print(logparams);
	vx1:=LinearAlgebra[DeleteColumn](logtrace,[2,3,4,5]);
	vx2:=LinearAlgebra[DeleteColumn](logtrace,[1,3,4,5]);
	vx3:=LinearAlgebra[DeleteColumn](logtrace,[1,2,4,5]);
	vx4:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,5]);
	vx5:=LinearAlgebra[DeleteColumn](logtrace,[1,2,3,4]);
	pol1:=logparams[1]+~vx3;
	#k1x3
	#print(pol1);
	pol2:=logparams[2]+~vx1;
	#k2x1
	#print(pol2);
	pol3:=logparams[3]+~vx2;
	#k3x2
	pol4:=logparams[4]+~vx2+~vx5;
	#k4x2x5
	pol5:=map(x->logparams[5],vx1);
	#k6
	#print(pol5);
	pol6:=logparams[6]+~vx3+~vx3+~vx4;
	#k9x3^2x4
	pol7:=logparams[7]+~vx4;
	#k10x4
	n:=LinearAlgebra[RowDimension](d);
	#orders:=[];
	listeq:=[];
	eq1:="";
	cpt:=0;
	for i from 1 to n do
		orders:=[min(pol1[i][1],pol2[i][1],pol3[i][1])-vx1[i][1], min(pol2[i][1],pol3[i][1],pol4[i][1])-vx2[i][1], min(pol1[i][1],pol6[i][1],pol7[i][1])-vx3[i][1], min(pol4[i][1],pol6[i][1],pol7[i][1])-vx4[i][1], min(pol5[i][1], pol4[i][1])-vx5[i][1]];
		eq2:=CheckEquilibrationAndBranchesd1e11([round(prec*pol1[i][1])/prec,round(prec*pol2[i][1])/prec,round(prec*pol3[i][1])/prec,round(prec*pol4[i][1])/prec,round(prec*pol5[i][1])/prec,round(prec*pol6[i][1])/prec,round(prec*pol7[i][1])/prec],orders,[vx1[i][1],vx2[i][1],vx3[i][1],vx4[i][1],vx5[i][1]]);
		#print(eq2);
		if (eq1<>eq2) then
			listeq:=[op(listeq),[eq1,cpt]];
			eq1:=eq2;
			cpt:=1;
		else
			cpt:=cpt+1;
		end if;
	od;
	listeq:=[op(listeq),[eq1,cpt]];
	return listeq;
end:

TropicalTransitions:=proc(eps,prec,params)
	#This is a work function, not a general one. Return the list of transitions we have between partial equilibrations.
	local i, j, filelist, transitions, trace, stringfile, trans;
	filelist:=[];
	transitions:={};

	for i from 1 to 9 do
    	stringfile:=cat("./00",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 10 to 99 do
    	stringfile:=cat("./0",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 100 to 375 do
    	stringfile:=cat("./",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 1 to 375 do
	    trace:=TraceAnalyserd1e11(filelist[i],eps,prec,params);
	    #print(trace);
	    for j from 1 to numelems(trace)-1 do
	    	trans:={[trace[j][1],trace[j+1][1]]};
	    	transitions:=transitions union trans;
	    od;
	od;
	return transitions;
end:

TropicalCountedTransitions:=proc(eps,prec,params)
	#This is a work function, not a general one. Return the list of transitions we have between partial equilibrations and how many time it happens.
	local i, j, k, notintransition, filelist, transitions, trace, stringfile, trans, nbtrans, diftrans;
	filelist:=[];
	transitions:=Array([]);
	notintransition:=true;
	nbtrans:=0;
	diftrans:=1;
	for i from 1 to 9 do
    	stringfile:=cat("./00",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 10 to 99 do
    	stringfile:=cat("./0",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 100 to 375 do
    	stringfile:=cat("./",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 1 to 375 do
	    trace:=TraceAnalyserd1e11(filelist[i],eps,prec,params);
	    #print(trace);
	    for j from 1 to numelems(trace)-1 do
	    	trans:=[trace[j][1],trace[j+1][1]];
	    	notintransition:=true;
	    	for k from 1 to numelems(transitions) do
	    		if (trans = transitions[k][1]) then
	    			transitions[k][2]:=transitions[k][2]+1;
	    			nbtrans:=nbtrans+1;
	    			notintransition:=false;
	    			break;
	    		end if;
	    	od;
	    	if notintransition then
	    		transitions(diftrans):=[trans,1];
	    		diftrans:=diftrans+1;
	    		nbtrans:=nbtrans+1;
	    	end if;
	    od;
	od;
	return [transitions,nbtrans];
end:

TropicalTraceTransitions:=proc(eps,prec,params)
	#This is a work function, not a general one. Return the list of transitions we have between partial equilibrations and how many time it happens (max 1 for a trace).
	local i, j, k, t, notintransition, filelist, transitions, trace, stringfile, diftrans, tracetransi, tracetransitions;
	filelist:=[];
	transitions:=Array([]);
	tracetransitions:={};
	notintransition:=true;
	diftrans:=1;
	for i from 1 to 9 do
    	stringfile:=cat("./00",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 10 to 99 do
    	stringfile:=cat("./0",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 100 to 375 do
    	stringfile:=cat("./",convert(i,string),".txt");
    	filelist:=[op(filelist),stringfile];
    od;
    for i from 1 to 375 do
    	tracetransitions:={};
	    trace:=TraceAnalyserd1e11(filelist[i],eps,prec,params);
	    #print(trace);
	    for j from 1 to numelems(trace)-1 do
	    	tracetransi:={[trace[j][1],trace[j+1][1]]};
	    	tracetransitions:=tracetransitions union tracetransi;
	    od;
    	#trans:=[trace[j][1],trace[j+1][1]];
    	for t in tracetransitions do
    		notintransition:=true;
	    	for k from 1 to numelems(transitions) do
	    		if (t = transitions[k][1]) then
	    			transitions[k][2]:=transitions[k][2]+1;
	    			notintransition:=false;
	    			break;
	    		end if;
	    	od;
    		if notintransition then
	    		transitions(diftrans):=[t,1];
	    		diftrans:=diftrans+1;
	    	end if;
	    od;
	od;
	return transitions;
end:

WriteBranchesInFile:=proc(file1,file2,transitions)
	#This is a work function, not a general one. Add the code nnnnnebb obtained by TraceAnalyserd1e11 to the given file as a new column.
	local linein, lineout, i, cnt, s;
	fclose(file1);
	fclose(file2);
	#fopen(file1);
	#fopen(file2);
	s:=0;
	#linein:=readline(file1);
	#lineout:=cat(linein," ","eqbranch");
	#s:=s+writeline(file2,convert(lineout,string));
	#linein:=readline(file1);
	#lineout:=cat(linein," ","na");
	#s:=s+writeline(file2,convert(lineout,string));
	cnt:=0;
	i:=1;
	while (linein <> 0) do
		while (i <= numelems(transitions)) do
			cnt:=transitions[i][2];
			while (cnt > 0) do
				linein:=readline(file1);
				lineout:=cat(linein," ",transitions[i][1]);
				s:=s+writeline(file2,convert(lineout,string));
				cnt:=cnt-1;
			od;
			i:=i+1;
		od;
		linein:=readline(file1);
		#lineout:=cat(linein," ",transitions[i][1]);
		#s:=s+writeline(file2,convert(lineout,string));
	od;
	fclose(file1);
	fclose(file2);
	return s;
end:

WriteMultipleBranchesInFiles:=proc(eps,prec,params)
	#This is a work function, not a general one. Add the code nnnnnebb obtained by TraceAnalyserd1e11 to the files as a new column.
	local i, filelist, trace, stringfile, stringfile2, w, filelist2;
	filelist:=[];
	filelist2:=[];
	for i from 1 to 9 do
    	stringfile:=cat("./00",convert(i,string),"n.txt");
    	filelist:=[op(filelist),stringfile];
    	stringfile2:=cat("./00",convert(i,string),"nbis.txt");
    	filelist2:=[op(filelist2),stringfile2];
    od;
    for i from 10 to 99 do
    	stringfile:=cat("./0",convert(i,string),"n.txt");
    	filelist:=[op(filelist),stringfile];
    	stringfile2:=cat("./0",convert(i,string),"nbis.txt");
    	filelist2:=[op(filelist2),stringfile2];
    od;
    for i from 100 to 375 do
    	stringfile:=cat("./",convert(i,string),"n.txt");
    	filelist:=[op(filelist),stringfile];
    	stringfile2:=cat("./",convert(i,string),"nbis.txt");
    	filelist2:=[op(filelist2),stringfile2];
    od;
    for i from 1 to 375 do
    	print(i);
	    trace:=TraceAnalyserd1e11(filelist[i],eps,prec,params);
	    w:=WriteBranchesInFile(filelist[i],filelist2[i],trace);
	od;
	return w;
end:

#Do not forget to use with(PolyhedralSets) before read("<PATH>/Showfiles.mpl").
