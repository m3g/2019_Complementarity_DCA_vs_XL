# Requires:
#path("/home/leandro/programs/JPleXL")
#using JPleXL # from https://github.com/mcubeg/JPleXL
#using Printf

d(x1,x2,x3,y1,y2,y3) = sqrt( (x1-y1)^2 + (x2-y2)^2 + (x3-y3)^2 )

function radius(name)
  if name[1:1] == "H"
    return 1.00
  elseif name[1:1] == "C"
    return 2.00
  elseif name[1:1] == "N"
    return 1.85
  elseif name[1:1] == "O"
    return 1.70
  elseif name[1:1] == "S"
    return 2.00
  elseif name[1:1] == "P"
    return 2.15
  elseif name[1:2] == "MG"
    return 1.185
  elseif name[1:2] == "MG"
    return 1.185
  elseif name[1:2] == "FE"
    return 1.90
  else
    println(" Could not find radius for ", name, " using 1.5A ")
    return 1.5
  end
end

function surface_dist(pdb)

  meshsize = 0.5
  probe = 1.4
  nprobe = round(Int64,(probe+3.00)/meshsize)+1 

  xmin = zeros(3)
  xmax = zeros(3)
  xmin[1] = minimum( atom -> atom.x, pdb ) - probe
  xmin[2] = minimum( atom -> atom.y, pdb ) - probe
  xmin[3] = minimum( atom -> atom.z, pdb ) - probe
  xmax[1] = maximum( atom -> atom.x, pdb ) + probe
  xmax[2] = maximum( atom -> atom.y, pdb ) + probe
  xmax[3] = maximum( atom -> atom.z, pdb ) + probe

  ngrid = zeros(Int64,3)
  @. ngrid = round((xmax - xmin) / meshsize)
  println(ngrid)

  xsize = zeros(3)
  @. xsize = (xmax-xmin)/ngrid
  println(xsize)

  grid = Array{Float64}(undef,ngrid[1],ngrid[2],ngrid[3],3)
  hasatom = Array{Bool}(undef,ngrid[1],ngrid[2],ngrid[3]) 
  for i in 1:ngrid[1]
    for j in 1:ngrid[2]
       for k in 1:ngrid[3]
         grid[i,j,k,1] = xmin[1] + (i-1)*xsize[1]
         grid[i,j,k,2] = xmin[2] + (j-1)*xsize[2]
         grid[i,j,k,3] = xmin[3] + (k-1)*xsize[3]
         hasatom[i,j,k] = false
       end
    end
  end
  for atom in pdb
    ix = round(Int64,(atom.x - xmin[1])/xsize[1])
    iy = round(Int64,(atom.y - xmin[2])/xsize[2])
    iz = round(Int64,(atom.z - xmin[3])/xsize[3])
    for i in max(1,ix-nprobe):min(ix+nprobe,ngrid[1])
      for j in max(1,iy-nprobe):min(iy+nprobe,ngrid[2])
        for k in max(1,iz-nprobe):min(iz+nprobe,ngrid[3])
          dist = d(atom.x,atom.y,atom.z,grid[i,j,k,1],grid[i,j,k,2],grid[i,j,k,3])
          if dist < probe + radius(atom.name)
            hasatom[i,j,k] = true
          end
        end
      end
    end
  end

  dmin = Vector{Float64}(undef,length(pdb))
  @. dmin = 10000.

  iatom = 0
  for atom in pdb
    iatom = iatom + 1
    ix = round(Int64,(atom.x - xmin[1])/xsize[1])
    iy = round(Int64,(atom.y - xmin[2])/xsize[2])
    iz = round(Int64,(atom.z - xmin[3])/xsize[3])
    found = false
    np = nprobe + 5
    while ! found
      for i in max(1,ix-np):min(ix+np,ngrid[1])
        for j in max(1,iy-np):min(iy+np,ngrid[2])
          for k in max(1,iz-np):min(iz+np,ngrid[3])
            if ! hasatom[i,j,k]
              dist = d(atom.x,atom.y,atom.z,grid[i,j,k,1],grid[i,j,k,2],grid[i,j,k,3])
              dmin[iatom] = min(dmin[iatom],dist-radius(atom.name))
              found = true
            end
          end
        end
      end
      np = np + 10
    end
  end
 
  ires = 0
  nres = 0
  for atom in pdb  
    if atom.resnum != ires
      nres = nres + 1
      ires = atom.resnum
    end
  end
  println(nres)

  dmin_res = Vector{Float64}(undef,nres) 
  ires = 0
  nres = 0
  iatom = 0
  for atom in pdb
    iatom = iatom + 1
    if atom.resnum != ires
      ires = atom.resnum
      nres = nres + 1
      dmin_res[nres] = dmin[iatom]
    else
      dmin_res[nres] = min(dmin[iatom],dmin_res[nres])
    end
  end

  dmin_res_atom = similar(dmin)
  iatom = 0
  ires = 0
  nres = 0
  for atom in pdb
    iatom = iatom + 1
    if atom.resnum != ires
      ires = atom.resnum
      nres = nres + 1
    end
    dmin_res_atom[iatom] = dmin_res[nres]
  end

  return dmin_res, dmin, dmin_res_atom

end

# Example:
#pdb = JPleXL.readPDB("../../../Correlation/Data/1AMM_clean.pdb")
#dmin_res, dmin, dmin_atom = surface_dist(pdb)


