
dcac = "red"
xlc = "blue"

plot(layout=(2,2))

# Mutual information plots
sp = 1

type = "cc"  
type = "p"

if type == "cc"
  plot!(xlabel=L"\textrm{\sffamily Contact Distance Threshold / \AA}")
  plot!(ylabel=L"\textrm{\sffamily Mutual Overlap}")
  for i in 1:length(structs)
    c = get(ColorSchemes.rainbow,i./length(structs))
    plot!(subplot=1,dcontact,all_cc_rand[i,:],linewidth=2,label=structs[i],color=c)
    plot!(subplot=1,dcontact,[ all_cc_obs[i] for j in 1:ndcontact ],label="",linewidth=2,color=c,linestyle=:dot)
  end
  annotate!([( 9.1,  0.03, text("DCA x XL",12,:center) )],subplot=1)
  annotate!([( 9.5,  0.60, text("Random",12,:center) )],subplot=1)
end
if type == "p"
  plot!(xlabel=L"\textrm{\sffamily Contact Distance Threshold / \AA}")
  plot!(ylabel=L"\textrm{\sffamily Probability of~}n\leq n_{XL\cap DCA}")
  for i in 1:length(structs)
    c = get(ColorSchemes.rainbow,i./length(structs))
    plot!(subplot=1,dcontact,all_p_rand[i,:],linewidth=2,label=structs[i],color=c)
  end
  #annotate!([( 9.1,  0.03, text("DCA x XL",12,:center) )],subplot=1)
  #annotate!([( 9.5,  0.60, text("Random",12,:center) )],subplot=1)
end

#order = sort([ i for i in 1:length(structs) ], by = i -> sizes[i])
#all_cc_obs = all_cc_obs[order]
#all_cc_rand = all_cc_rand[order]
#sort!(sizes)
#histogram!(all_cc_obs,label="observed MO",nbins=6,alpha=0.5,color="yellow",subplot=sp)
#histogram!(all_cc_rand,label="random MO",nbins=12,alpha=0.5,color="darkgreen",subplot=sp)
#bar!(xlabel=L"\textrm{\sffamily Mutual Overlap}",subplot=sp)
#bar!(ylabel=L"\textrm{\sffamily Count}",subplot=sp)

annotate!([( -1.8-16.5,  500, "A" )],fontsize=48,subplot=4)
annotate!([( -1.8,  500, "B" )],fontsize=48,subplot=4)
annotate!([( -1.8-16.5, 200, "C" )],fontsize=48,subplot=4)
annotate!([( -1.8, 200, "D" )],fontsize=48,subplot=4)

# residue type plot
sp = 2

bar!(all_types_dca,alpha=0.5,xrotation=60,label="DCAs",color=dcac,subplot=sp)
bar!(all_types_xl,alpha=0.5,xrotation=60,label="XLs",xticks=(1:1:20,restypes),color=xlc,subplot=sp)
bar!(xlabel=L"\textrm{\sffamily Residue Type}",ylabel=L"\textrm{\sffamily Count}",subplot=sp)

# CA - CA distances 
sp = 3

xl_bin = ( maximum(all_xl) - minimum(all_xl) ) / 40
ndcabins = round(Int64,( maximum(all_dca) - minimum(all_dca) ) / xl_bin)
histogram!(all_dca,bins=ndcabins,label="DCAs",alpha=0.5,color=dcac,subplot=sp)
histogram!(all_xl,bins=40,label="XLs",alpha=1.0,color=xlc,subplot=sp)
histogram!(xlabel=L"\textrm{\sffamily C}\alpha\textrm{\sffamily~Euclidean Distance} / \textrm{\sffamily~\AA}",subplot=sp)
histogram!(ylabel=L"\textrm{\sffamily Count}",subplot=sp)
m1 = mean(all_dca)
m2 = mean(all_xl)
scatter!([m1,m1,m1,m1],[100,104,108,112],label="",color=dcac,linewidth=3,linestyle=:dot,subplot=sp,markersize=3)
scatter!([m2,m2,m2,m2],[100,104,108,112],label="",color=xlc,linewidth=3,linestyle=:dot,subplot=sp,markersize=3)

# Maximum distance to surface
sp = 4

xl_bin = ( maximum(all_max_xl_surfdist) - minimum(all_max_xl_surfdist) ) / 40
ndcabins = round(Int64,( maximum(all_max_dca_surfdist) - minimum(all_max_dca_surfdist) ) / xl_bin)
histogram!(all_max_dca_surfdist,bins=ndcabins,label="DCAs",alpha=0.5,color=dcac,subplot=sp)
histogram!(all_max_xl_surfdist,bins=40,label="XLs",alpha=1.0,color=xlc,subplot=sp)
histogram!(xlabel=L"\textrm{\sffamily Maximum C}\alpha\textrm{\sffamily~distance to surface /}\textrm{\sffamily~\AA}",subplot=sp)
histogram!(ylabel=L"\textrm{\sffamily Count}",subplot=sp)

x, y = M3GTools.density(all_max_contact_surfdist,step=1.0,vmin=1.0)
scale = 205/maximum(y)
y = scale*y
plot!(x,y,subplot=sp,label="All contacts",linewidth=2,color="green",alpha=0.8)

m1 = mean(all_max_dca_surfdist)
m2 = mean(all_max_xl_surfdist)
scatter!([m1,m1,m1,m1],[180,188,196,204],label="",color=dcac,linewidth=3,linestyle=:dot,subplot=sp,markersize=3)
scatter!([m2,m2,m2,m2],[180,188,196,204],label="",color=xlc,linewidth=3,linestyle=:dot,subplot=sp,markersize=3)

plot!(size=(750,750))
savefig("./all.pdf")
#savefig("./all.png")

# Distances to surface
plot()

xl_bin = ( maximum(all_xl_surfdist) - minimum(all_xl_surfdist) ) / 40
ndcabins = round(Int64,( maximum(all_dca_surfdist) - minimum(all_dca_surfdist) ) / xl_bin)
histogram!(all_dca_surfdist,bins=ndcabins,label="DCAs",alpha=0.5,color=dcac)
histogram!(all_xl_surfdist,bins=40,label="XLs",alpha=1.0,color=xlc)
histogram!(xlabel=L"\textrm{\sffamily C}\alpha\textrm{\sffamily~distance to surface /}\textrm{\sffamily~\AA}")
histogram!(ylabel=L"\textrm{\sffamily Count}")

#histogram!(all_contact_surfdist,bins=40,label="XLs",alpha=1.0,color="green")

x, y = M3GTools.density(all_contact_surfdist,vmin=1,vmax=12.5,step=0.5,nbins=25)
scale = 400/maximum(y)
y = scale*y
plot!(x,y,label="All contacts",linewidth=3,color="green",alpha=1.0)

m1 = mean(all_dca_surfdist)
m2 = mean(all_xl_surfdist)
m3 = mean(all_contact_surfdist)
scatter!([m1,m1,m1,m1],[180,188,196,204],label="",color=dcac,linewidth=3,linestyle=:dot,markersize=3)
scatter!([m2,m2,m2,m2],[180,188,196,204],label="",color=xlc,linewidth=3,linestyle=:dot,markersize=3)
scatter!([m3,m3,m3,m3],[180,188,196,204],label="",color="green",linewidth=3,linestyle=:dot,markersize=3)

plot!(size=(350,350))
savefig("./CA_to_surface.pdf")





