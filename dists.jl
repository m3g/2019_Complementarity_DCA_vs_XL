path("/home/leandro/programs/M3GTools")
path("/home/leandro/programs/JPleXL")
nogtk()
using JPleXL
using M3GTools
using Plots
using LaTeXStrings
using Printf
using Statistics
using ColorSchemes
ENV["GKS_ENCODING"]="utf-8"
#pyplot()

include("./surface.jl")

restypes = [ "ALA",
             "ARG",
             "ASN",
             "ASP",
             "CYS",
             "GLU",
             "GLN",
             "GLY",
             "HIS",
             "ILE",
             "LEU",
             "LYS",
             "MET",
             "PHE",
             "PRO",
             "SER",
             "THR",
             "TRP",
             "TYR",
             "VAL" ]


basedir = "/home/leandro/Drive/Alunos/Ricardo/Proteins/review/CompleteDCA_Ricardo"
list_full_DCA_Ricardo = [ 
                          [ 162,   basedir*"/1AMM/1AMM_PF00030_align_cleaned_ranked_matched.DI" ],
                          [ 194,   basedir*"/1ARB/1ARB_PF00089_align_cleaned_ranked_matched.DI" ],
                          [ 226,   basedir*"/1ATG/1ATG_PF13531_align_cleaned_ranked_matched.DI" ],
                          [ 106,   basedir*"/1B0B/1B0B_PF00042_align_cleaned_ranked_matched.DI" ],
                          [ 306,   basedir*"/1BXO/1BXO_PF00026_aligned_cleaned_ranked_matched.DI" ],
                          [ 89 ,   basedir*"/1C52/1C52_PF00034_align_ranked_matched.DI" ],
                          [ 66 ,   basedir*"/1C75/1C75_PF13442_align_cleaned_ranked_matched.DI" ],
                          [ 112,   basedir*"/1D06/1D06_PF00989_align_cleaned_ranked_matched.DI" ],
                          [ 80 ,   basedir*"/1D4T/1D4T_PF00017_align_cleaned_ranked_matched.DI" ],
                          [ 102,   basedir*"/1EW4/1EW4_PF01491_align_ranked_matched.DI" ],
                          [ 86 ,   basedir*"/1FK5/1FK5_PF00234_align_cleaned_ranked_matched.DI" ],
                          [ 184,   basedir*"/1G67/1G67_PF02581_align_cleaned_ranked_matched.DI" ],
                          [ 52 ,   basedir*"/1G6X/1G6X_PF00014_align_cleaned_ranked_matched.DI" ],
                          [ 221,   basedir*"/PDBs/1G8A/1G8A_PF01269_align_cleaned_ranked_matched.DI" ],
                          [ 232,   basedir*"/1GCI/1GCI_PF00082_align_cleaned_ranked_matched.DI" ] 
                        ]


dcontact = [ 8., 9., 10., 11., 12., 13., 14., 15. ]
ndcontact = length(dcontact)

function inpdb(i,j,pdb)
  if findfirst( atom -> atom.name == "CA" && atom.resnum == i , pdb ) == nothing 
    return false
  end
  if findfirst( atom -> atom.name == "CA" && atom.resnum == j , pdb ) == nothing 
    return false
  end
  return true
end

function addtype!(name,types)
  try 
    ires = findfirst( x -> x == name, restypes )
    types[ires] = types[ires] + 1 
  catch
    println(name)
  end
end

function repeated(x)
  for i in 1:length(x)
    if count( a -> a == x[i], x ) > 1
      return i
    end
  end
  return 0
end

function rand_no_repeat(N,n)
  x = rand(1:N,n)
  repeat = repeated(x)
  while repeat != 0
    x[repeat] = rand(1:N)
    repeat = repeated(x)
  end
  return x
end

function cc(nca,nxl,nequal)
  nmin = min(nxl,nca)
  #cc = ( nequal - ( nmin - nequal ) ) / nmin 
  cc = nequal / nmin 
  return cc
end

d(x,y) = sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3]-y[3])^2 )
  
function number_of_contacts(cas,ncas,dcon)
  nc = 0
  for i in 1:ncas-5
    for j in i+5:ncas
      dij = d(cas[i,:],cas[j,:]) 
      if dij < dcon
        nc = nc + 1
      end
    end
  end
  return nc
end

function distances(name,linkdata,pdb)

  # Read cross-linking data

  file = open(linkdata,"r")
  nxl = 0 
  nxl_all = 0 
  for line in eachline(file)
    data = split(line)
    if ( length(data) == 4 )
      i = parse(Int64,data[1])
      j = parse(Int64,data[2])
      if inpdb(i,j,pdb)
        nxl_all = nxl_all + 1
        if abs(i-j) > 4
          nxl = nxl + 1
        end
      else
        println(" Residues not found: ", i," ", j)
      end
    end
  end
  close(file)
  println(" nxl with i-j < 4 = ", nxl_all )
  println(" nxl = ", nxl )

  xl = Matrix{Int64}(undef,nxl,2)
  file = open(linkdata,"r")
  ixl = 0
  for line in eachline(file)
    data = split(line)
    if ( length(data) == 4 )
      i = parse(Int64,data[1])
      j = parse(Int64,data[2])
      if inpdb(i,j,pdb)
        if abs(i-j) > 4
          ixl = ixl + 1
          xl[ixl,1] = i
          xl[ixl,2] = j
        end
      end
    end
  end
  close(file)

  # Reading full DCA data from list_full_DCA_Ricardo 

  dcafile = Nothing
  ndca = 0
  for file in list_full_DCA_Ricardo
     if occursin(name,file[2])
      ndca = file[1]
      dcafile = file[2]
    end
  end
  println(dcafile)
  println(" ndca = ", ndca )
  dca = Matrix{Int64}(undef,ndca,2)
  ndca_plus = nxl
  dca_plus = Matrix{Int64}(undef,ndca_plus,2)
  file = open(dcafile,"r")
  idca = 0
  for line in eachline(file)
    data = split(line)
    if length(data) == 2 
      i = parse(Int64,data[1])
      j = parse(Int64,data[2])
      if inpdb(i,j,pdb)
        if abs(i-j) > 4 
          idca = idca + 1
          if idca <= ndca
            dca[idca,1] = i
            dca[idca,2] = j
          elseif idca <= ndca + nxl
            dca_plus[idca-ndca,1] = i
            dca_plus[idca-ndca,2] = j
          else
            break
          end
        end
      else
        println(" Residues not found, DCA: ", i," ",j)
      end 
    end
  end
  close(file)

  # Count number of common contacts

  nequal = 0
  for i in 1:nxl
    for j in 1:ndca
      if ( xl[i,1] == dca[j,1] && xl[i,2] == dca[j,2] ) ||
         ( xl[i,1] == dca[j,2] && xl[i,2] == dca[j,1] ) 
        nequal = nequal + 1
      end
    end
  end

  xca_xl = [ [ zeros(3) for i in 1:2 ] for j in 1:nxl ]
  xca_dca = [ [ zeros(3) for i in 1:2 ] for j in 1:ndca ]
  xca_dca_plus = [ [ zeros(3) for i in 1:2 ] for j in 1:ndca_plus ]

  types_xl = zeros(Int64,20)
  types_dca = zeros(Int64,20)
  
  # Number of CAs
  cas = JPleXL.xCA(pdb)
  ncas = size(cas)[1]

  dsurf_res, dsurf_atom, dsurf_res_atom = surface_dist(pdb)

  dca_surfdist = Vector{Float64}(undef,ndca)
  xl_surfdist = Vector{Float64}(undef,nxl)
  max_xl_surfdist = zeros(Float64,nxl)
  max_dca_surfdist = zeros(Float64,ndca)
  all_ca_surfdist = Vector{Float64}(undef,ncas)

  found_xl = zeros(Bool,nxl,2)
  found_dca = zeros(Bool,ndca,2)

  nfound_xl = 0
  nfound_dca = 0
  av_CA_surfdist = 0.
  iCA = 0
  for ipdb in 1:length(pdb)
    if pdb[ipdb].name == "CA"
      av_CA_surfdist = av_CA_surfdist + dsurf_atom[ipdb] 
      iCA = iCA + 1
      all_ca_surfdist[iCA] = dsurf_atom[ipdb]
      for i in 1:nxl
        if pdb[ipdb].resnum == xl[i,1]
          xca_xl[i][1][1] = pdb[ipdb].x  
          xca_xl[i][1][2] = pdb[ipdb].y  
          xca_xl[i][1][3] = pdb[ipdb].z  
          nfound_xl = nfound_xl + 1
          addtype!(pdb[ipdb].resname,types_xl)
          found_xl[i,1] = true
          xl_surfdist[i] = dsurf_atom[ipdb]
          max_xl_surfdist[i] = max(max_xl_surfdist[i],dsurf_atom[ipdb])  
        end
        if pdb[ipdb].resnum == xl[i,2]
          xca_xl[i][2][1] = pdb[ipdb].x  
          xca_xl[i][2][2] = pdb[ipdb].y  
          xca_xl[i][2][3] = pdb[ipdb].z  
          nfound_xl = nfound_xl + 1
          addtype!(pdb[ipdb].resname,types_xl)
          found_xl[i,2] = true
          xl_surfdist[i] = dsurf_atom[ipdb]
          max_xl_surfdist[i] = max(max_xl_surfdist[i],dsurf_atom[ipdb])  
        end
      end
      for i in 1:ndca
        if pdb[ipdb].resnum == dca[i,1]
          xca_dca[i][1][1] = pdb[ipdb].x  
          xca_dca[i][1][2] = pdb[ipdb].y  
          xca_dca[i][1][3] = pdb[ipdb].z  
          nfound_dca = nfound_dca + 1
          addtype!(pdb[ipdb].resname,types_dca)
          found_dca[i,1] = true
          dca_surfdist[i] = dsurf_atom[ipdb]
          max_dca_surfdist[i] = max(max_dca_surfdist[i],dsurf_atom[ipdb])  
        end
        if pdb[ipdb].resnum == dca[i,2]
          xca_dca[i][2][1] = pdb[ipdb].x  
          xca_dca[i][2][2] = pdb[ipdb].y  
          xca_dca[i][2][3] = pdb[ipdb].z  
          nfound_dca = nfound_dca + 1
          addtype!(pdb[ipdb].resname,types_dca)
          found_dca[i,2] = true
          dca_surfdist[i] = dsurf_atom[ipdb]
          max_dca_surfdist[i] = max(max_dca_surfdist[i],dsurf_atom[ipdb])  
        end
      end
      for i in 1:ndca_plus
        if pdb[ipdb].resnum == dca_plus[i,1]
          xca_dca_plus[i][1][1] = pdb[ipdb].x  
          xca_dca_plus[i][1][2] = pdb[ipdb].y  
          xca_dca_plus[i][1][3] = pdb[ipdb].z  
        end
        if pdb[ipdb].resnum == dca_plus[i,2]
          xca_dca_plus[i][2][1] = pdb[ipdb].x  
          xca_dca_plus[i][2][2] = pdb[ipdb].y  
          xca_dca_plus[i][2][3] = pdb[ipdb].z  
        end
      end
    end
  end
  av_CA_surfdist = av_CA_surfdist / ncas
  println(" Aveage CA distance to surface: ", av_CA_surfdist, " length = ", ncas)
  for i in 1:nxl
    if ! found_xl[i,1] 
      println(" XL not found: ", xl[i,1] ) 
    end
    if ! found_xl[i,2]
      println(" XL not found: ", xl[i,2] ) 
    end
  end
  for i in 1:ndca
    if ! found_dca[i,1] 
      println(" DCA not found: ", dca[i,1] ) 
    end
    if ! found_dca[i,2]
      println(" DCA not found: ", dca[i,2] ) 
    end
  end
  
  d_xl = Vector{Float64}(undef,nxl)
  d_dca = Vector{Float64}(undef,ndca)
  d_dca_plus = Vector{Float64}(undef,ndca_plus)
  
  for i in 1:nxl
    x = xca_xl[i][1]
    y = xca_xl[i][2]
    d_xl[i] = d(x,y)
  end
  
  for i in 1:ndca_plus
    x = xca_dca_plus[i][1]
    y = xca_dca_plus[i][2]
    d_dca_plus[i] = d(x,y)
  end

  tol = 8.
  n_true_dca = 0
  n_true_dca_12 = 0
  for i in 1:ndca
    x = xca_dca[i][1]
    y = xca_dca[i][2]
    d_dca[i] = d(x,y)
    if ( d_dca[i] < tol )
      n_true_dca = n_true_dca + 1
    end
    if ( d_dca[i] < 12. )
      n_true_dca_12 = n_true_dca_12 + 1
    end
  end

  # Print XLs to file
  file = open("./constraints_final/"*name*"_XL.dat","w")
  for i in 1:nxl
    println(file,@sprintf("%8i %8i %10.3f",xl[i,1],xl[i,2],d_xl[i]))
  end
  close(file)
  # Print DCAs to file
  file = open("./constraints_final/"*name*"_DCA.dat","w")
  for i in 1:ndca
    println(file,@sprintf("%8i %8i %10.3f",dca[i,1],dca[i,2],d_dca[i]))
  end
  close(file)
  # Print DCA_plus to file
  file = open("./constraints_final/"*name*"_DCA_plus.dat","w")
  for i in 1:ndca_plus
    println(file,@sprintf("%8i %8i %10.3f",dca_plus[i,1],dca_plus[i,2],d_dca_plus[i]))
  end
  close(file)

  
  cc_obs = cc(ndca,nxl,nequal)
  cc_rand = Vector{Float64}(undef,ndcontact)
  av_nequal = Vector{Float64}(undef,ndcontact)
  p_random = Vector{Float64}(undef,ndcontact)
  n_random = Vector{Int64}(undef,ndcontact)
  f_size = Vector{Float64}(undef,ndcontact)

  for id in 1:ndcontact
    dcon = dcontact[id]

    # Number of contacts 
    nc = number_of_contacts(cas,ncas,dcon)
    f_size[id] = nc / ncas

    # Computing random result
    n_random[id] = 0
    av_nequal[id] = 0.
    ntrial = 1000
    for i in 1:ntrial
      random_dca = rand_no_repeat(nc,ndca)
      random_xl = rand_no_repeat(nc,nxl)
      nequal_random = 0
      for i in 1:nxl
        nequal_random = nequal_random + count( x -> x == random_xl[i], random_dca )
      end
      if nequal_random <= nequal 
        n_random[id] = n_random[id] + 1
      end
      av_nequal[id] = av_nequal[id] + nequal_random
    end
    av_nequal[id] = av_nequal[id] / ntrial
    p_random[id] = n_random[id] / ntrial
    cc_rand[id] = cc(ndca,nxl,av_nequal[id])

  end

  # All-contact distances to surface
  #nc_all = Int((ncas*ncas-ncas)/2) - (ncas-1) - (ncas-2) - (ncas-3) - (ncas-4) 
  dcon = 12.
  nc = number_of_contacts(cas,ncas,dcon)
  max_contact_surfdist = Vector{Float64}(undef,nc)
  ic = 0
  for i in 1:ncas-5
    for j in i+5:ncas
      dij = d(cas[i,:],cas[j,:]) 
      if dij < dcon
        ic = ic + 1
        max_contact_surfdist[ic] = max(all_ca_surfdist[i],all_ca_surfdist[j])
      end
    end
  end

  println("-------------------------------------------------------------")
  println( name )
  println(" Number of residues: ", ncas)
  println(" Number of contacts (",dcon,"A): ",nc)
  println(" Number of DCA contacts: ", ndca)
  println(" N true DCA (", tol,"A): ", n_true_dca)
  println(" N true DCA (12A): ", n_true_dca_12)
  println(" Number of XL contacts: ", nxl)

  println(" Average nequal = ", av_nequal, " Observed = ", nequal )
  println(" Probability of random result = ", p_random*100, "%" )
  println(" cc_obs = ", cc_obs) 
  println(" cc_rand = ", cc_rand)
  println(" nc / ncas = ", nc / ncas )
  println("-------------------------------------------------------------")
 
  return d_xl, d_dca, types_xl, types_dca, 
         xl_surfdist, dca_surfdist, all_ca_surfdist, 
         cc_obs, cc_rand, p_random, f_size,
         max_xl_surfdist, max_dca_surfdist, max_contact_surfdist

end

#voltar
structs = [ "1AMM",
            "1ARB",
            "1ATG",
            "1B0B",
            "1BXO",
            "1C52",
            "1C75",
            "1D06",
            "1D4T",
            "1EW4",
            "1FK5",
            "1G67",
            "1G6X",
            "1G8A",
            "1GCI" ]

order = sort!([ i for i in 1:length(structs)],by=x->list_full_DCA_Ricardo[x][1])
structs = structs[order]

#structs = [ "1G8A", "1GCI" ]
#structs = [ "1C52" ]
#structs = [ "1D06" ]
#structs = [ "1GCI" ]
#structs = [ "1FK5" ]


function unionfile(s)
  dir = "../$s"
  files = readdir(dir)
  for file in files
    if occursin("UNION",file)
      return file
    end
  end
end

all_xl = Vector{Float64}(undef,0)
all_dca = Vector{Float64}(undef,0)
all_types_xl = zeros(Int64,20)
all_types_dca = zeros(Int64,20)
all_xl_surfdist = Vector{Float64}(undef,0)
all_dca_surfdist = Vector{Float64}(undef,0)
all_cc_obs = Vector{Float64}(undef,length(structs))
all_cc_rand = Matrix{Float64}(undef,length(structs),ndcontact)
all_p_rand = Matrix{Float64}(undef,length(structs),ndcontact)
all_f_size = Matrix{Float64}(undef,length(structs),ndcontact)
all_xl_surfdist = Vector{Float64}(undef,0)
all_dca_surfdist = Vector{Float64}(undef,0)
all_contact_surfdist = Vector{Float64}(undef,0)
all_max_xl_surfdist = Vector{Float64}(undef,0)
all_max_dca_surfdist = Vector{Float64}(undef,0)
all_max_contact_surfdist = Vector{Float64}(undef,0)
sizes = Vector{Float64}(undef,length(structs))
ipdb = 0
for s in structs
  println(s)
  global ipdb = ipdb + 1
  pdb = readPDB("../../../Correlation/Data/$s.pdb")
  sizes[ipdb] = length(filter( atom -> atom.name == "CA", pdb))
  uf = unionfile(s)
  linkdata = "../$s/$uf"
  d_xl, d_dca, types_xl, types_dca, 
       xl_surfdist, dca_surfdist, contact_surfdist,
       cc_obs, cc_rand, p_rand, f_size,
       max_xl_surfdist, max_dca_surfdist, max_contact_surfdist = distances(s,linkdata,pdb)
  all_cc_obs[ipdb] = cc_obs
  all_cc_rand[ipdb,:] = cc_rand
  all_p_rand[ipdb,:] = p_rand
  all_f_size[ipdb,:] = f_size
  append!(all_xl,d_xl)
  append!(all_dca,d_dca)
  append!(all_xl_surfdist,xl_surfdist)
  append!(all_dca_surfdist,dca_surfdist)
  append!(all_contact_surfdist,contact_surfdist)
  append!(all_max_xl_surfdist,max_xl_surfdist)
  append!(all_max_dca_surfdist,max_dca_surfdist)
  append!(all_max_contact_surfdist,max_contact_surfdist)
  @. all_types_xl = all_types_xl + types_xl
  @. all_types_dca = all_types_dca + types_dca
end

finalxl = open("xls.dat","w")
for i in 1:length(all_xl)
  println(finalxl,all_xl[i])
end
close(finalxl)
finaldca = open("dcas.dat","w")
for i in 1:length(all_dca)
  println(finaldca,all_dca[i])
end
close(finaldca)

include("./all_plots.jl")

