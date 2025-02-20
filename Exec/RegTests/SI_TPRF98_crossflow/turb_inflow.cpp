#include "turb_inflow.H"

void
read_turbulent_bc(int ix_limit[2], amrex::Real v_mean){
        amrex::Print() << "Reading turbulence from: " <<  PeleC::prob_parm_host->turb_ic.c_str() << std::endl;

        //get array sizes
        int z_coord = 1;
        get_inputs_turb(PeleC::prob_parm_host->turb_ic,z_coord);
        
        //define local variable containing size info
        int nx        = ix_limit[1] - ix_limit[0]; // number of points read from the file
        int nx_global = PeleC::prob_parm_host->nx_tmp; // this is the total number of points in x contained in the read file
        int ny        = PeleC::prob_parm_host->ny_tmp;
        int nz        = PeleC::prob_parm_host->nz_tmp;
        int nscal     = PeleC::prob_parm_host->nscal_tmp; 
        
        amrex::Print() << "Reading only: " << nscal << " scalars with nx =  " << nx 
                              << " and ny = " << ny << " and nz = " << nz << " points \n" <<std::endl;  

        //initialize variables for interpolation later on
        PeleC::h_prob_parm_device->ires_x = nx;
        PeleC::h_prob_parm_device->ires_y = ny;
        PeleC::h_prob_parm_device->ires_z = nz;

        amrex::Vector<double> data_tmp(nx * ny * nz * nscal);

        amrex::Vector<double> gridx_input(nx);
        amrex::Vector<double> gridy_input(ny);
        amrex::Vector<double> gridz_input(nz);

        // //get data from file
        read_turb_input_file(PeleC::prob_parm_host->turb_ic,
                        z_coord,
                        nx_global,
                        ny,
                        nz,
                        nscal,
                        ix_limit, // array containing the indexes for the points in x to be read
                        gridx_input, gridy_input, gridz_input, data_tmp);
        
        //copy data to the working arrays
        PeleC::prob_parm_host->T_xinput.resize(PeleC::h_prob_parm_device->ires_x);
        for (int i = 0; i < PeleC::h_prob_parm_device->ires_x; i++) {
          PeleC::prob_parm_host->T_xinput[i] = gridx_input[i]/v_mean; //transforming x-coord into time coordinate
        }
        PeleC::prob_parm_host->T_yinput.resize(PeleC::h_prob_parm_device->ires_y);
        for (int i = 0; i < PeleC::h_prob_parm_device->ires_y; i++) {
          PeleC::prob_parm_host->T_yinput[i] = gridy_input[i];
        }
        PeleC::prob_parm_host->T_zinput.resize(PeleC::h_prob_parm_device->ires_z);
        for (int i = 0; i < PeleC::h_prob_parm_device->ires_z; i++) {
          PeleC::prob_parm_host->T_zinput[i] = gridz_input[i];
        }
        PeleC::prob_parm_host->data_turb.resize(nx*ny*nz*nscal);
        for (int i = 0; i < nx*ny*nz*nscal; i++) {
          PeleC::prob_parm_host->data_turb[i] = data_tmp[i];
        }


        // Get the array differences.
        PeleC::prob_parm_host->T_xdiff.resize(PeleC::h_prob_parm_device->ires_x);
        std::adjacent_difference(
          PeleC::prob_parm_host->T_xinput.begin(), PeleC::prob_parm_host->T_xinput.end(),
          PeleC::prob_parm_host->T_xdiff.begin());
        PeleC::prob_parm_host->T_xdiff[0] = PeleC::prob_parm_host->T_xdiff[1];

        PeleC::prob_parm_host->T_ydiff.resize(PeleC::h_prob_parm_device->ires_y);
        std::adjacent_difference(
          PeleC::prob_parm_host->T_yinput.begin(), PeleC::prob_parm_host->T_yinput.end(),
          PeleC::prob_parm_host->T_ydiff.begin());
        PeleC::prob_parm_host->T_ydiff[0] = PeleC::prob_parm_host->T_ydiff[1];

        PeleC::prob_parm_host->T_zdiff.resize(PeleC::h_prob_parm_device->ires_z);
        std::adjacent_difference(
          PeleC::prob_parm_host->T_zinput.begin(), PeleC::prob_parm_host->T_zinput.end(),
          PeleC::prob_parm_host->T_zdiff.begin());
        PeleC::prob_parm_host->T_zdiff[0] = PeleC::prob_parm_host->T_zdiff[1];


        PeleC::prob_parm_host->data_ic.resize(PeleC::prob_parm_host->data.size());

        PeleC::prob_parm_host->xarray_turb.resize(PeleC::prob_parm_host->T_xinput.size());
        PeleC::prob_parm_host->xdiff_turb.resize(PeleC::prob_parm_host->T_xdiff.size());
        PeleC::prob_parm_host->yarray_turb.resize(PeleC::prob_parm_host->T_yinput.size());
        PeleC::prob_parm_host->ydiff_turb.resize(PeleC::prob_parm_host->T_ydiff.size());
        PeleC::prob_parm_host->zarray_turb.resize(PeleC::prob_parm_host->T_zinput.size());
        PeleC::prob_parm_host->zdiff_turb.resize(PeleC::prob_parm_host->T_zdiff.size());
        PeleC::prob_parm_host->data_ic_turb.resize(PeleC::prob_parm_host->data_turb.size());


        // Get pointer to the data
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_xinput.begin(),
        PeleC::prob_parm_host->T_xinput.end(),  PeleC::prob_parm_host->xarray_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_xdiff.begin(),
        PeleC::prob_parm_host->T_xdiff.end(), PeleC::prob_parm_host->xdiff_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_yinput.begin(),
        PeleC::prob_parm_host->T_yinput.end(),  PeleC::prob_parm_host->yarray_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_ydiff.begin(),
        PeleC::prob_parm_host->T_ydiff.end(), PeleC::prob_parm_host->ydiff_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_zinput.begin(),
        PeleC::prob_parm_host->T_zinput.end(),  PeleC::prob_parm_host->zarray_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->T_zdiff.begin(),
        PeleC::prob_parm_host->T_zdiff.end(), PeleC::prob_parm_host->zdiff_turb.begin()
          );

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->data_turb.begin(),
        PeleC::prob_parm_host->data_turb.end(), PeleC::prob_parm_host->data_ic_turb.begin()
          );


        PeleC::h_prob_parm_device->d_xarray_turb  = PeleC::prob_parm_host->xarray_turb.data();
        PeleC::h_prob_parm_device->d_xdiff_turb   = PeleC::prob_parm_host->xdiff_turb.data();
        PeleC::h_prob_parm_device->d_yarray_turb  = PeleC::prob_parm_host->yarray_turb.data();
        PeleC::h_prob_parm_device->d_ydiff_turb   = PeleC::prob_parm_host->ydiff_turb.data();
        PeleC::h_prob_parm_device->d_zarray_turb  = PeleC::prob_parm_host->zarray_turb.data();
        PeleC::h_prob_parm_device->d_zdiff_turb   = PeleC::prob_parm_host->zdiff_turb.data();
        PeleC::h_prob_parm_device->d_data_ic_turb = PeleC::prob_parm_host->data_ic_turb.data();


        // Dimensions of the input box
        PeleC::h_prob_parm_device->Lxturb =
          PeleC::prob_parm_host->T_xinput[nx - 1];// - PeleC::prob_parm_host->T_xinput[0];

        PeleC::h_prob_parm_device->Lyturb =
          PeleC::prob_parm_host->T_yinput[ny - 1] + 0.5 * PeleC::prob_parm_host->T_yinput[ny - 1];

        PeleC::h_prob_parm_device->Lzturb =
          PeleC::prob_parm_host->T_zinput[nz - 1] + 0.5 * PeleC::prob_parm_host->T_zinput[nz - 1];
}

void
get_inputs_turb(
  const std::string& iname,
  int z_coord)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile); //The problem is coming from here
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  std::string line;
  
  std::getline(iss, line);
  std::istringstream ninput(line);
  ninput >> PeleC::prob_parm_host->nscal_tmp;

  std::getline(iss, line);
  std::istringstream xinput(line);
  xinput >> PeleC::prob_parm_host->nx_tmp;
  
  std::getline(iss, line);
  std::istringstream yinput(line);
  yinput >> PeleC::prob_parm_host->ny_tmp;

  amrex::Print() << "In get_inputs_turb() reading data from: " <<  iname.c_str() << std::endl;
  if(z_coord){
    std::getline(iss, line);
    std::istringstream zinput(line);
    zinput >> PeleC::prob_parm_host->nz_tmp;
    PeleC::prob_parm_host->nscal_tmp = PeleC::prob_parm_host->nscal_tmp - 3; //remove x and y and z from count
    amrex::Print() << PeleC::prob_parm_host->nscal_tmp << " scalars found with nx =  " << PeleC::prob_parm_host->nx_tmp 
                          << " and ny = " << PeleC::prob_parm_host->ny_tmp << " and nz = " << PeleC::prob_parm_host->nz_tmp << " points \n" <<std::endl;
  }
  else{
    PeleC::prob_parm_host->nscal_tmp = PeleC::prob_parm_host->nscal_tmp - 2; //remove x and y from count
    PeleC::prob_parm_host->nz_tmp = 1;
    amrex::Print() << PeleC::prob_parm_host->nscal_tmp << " scalars found with nx =  " << PeleC::prob_parm_host->nx_tmp 
                          << " and ny = " << PeleC::prob_parm_host->ny_tmp << " points \n" << std::endl;
  }

};

void
read_turb_input_file(
  const std::string iname,
  int z_coord,
  int nx,
  int ny,
  int nz,
  int nscal,
  int ix_dir_limit[2],
  amrex::Vector<amrex::Real>& gridx_input,
  amrex::Vector<amrex::Real>& gridy_input,
  amrex::Vector<amrex::Real>& gridz_input,
  amrex::Vector<amrex::Real>& data)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  int i,j,k,n;

  amrex::Print() << "Reading new turbulence data with index: " << ix_dir_limit[0] << " " << ix_dir_limit[1] <<std::endl;

  std::istringstream iss(memfile);

  // Read the file
  int nlines = 0;
  std::string firstline, line;
  if(z_coord){
    for (i=0;i<nscal+7;i++){
      std::getline(iss, firstline); // skip header
    }  
  }
  else{
    for (i=0;i<nscal+5;i++){
      std::getline(iss, firstline); // skip header
    }    
  }


  // Read grid in X
  int cnt = 0;
  for (i=0;i<nx;i++){
    std::getline(iss, line);
    std::istringstream sinput(line);
    if(i >= ix_dir_limit[0] and i < ix_dir_limit[1]){
      sinput >> gridx_input[cnt];
      cnt++;
    }
  }

  // Read the data from the file
  cnt = 0;
  for (j=0;j<ny;j++){
    std::getline(iss, line);
    std::istringstream sinput(line);
    sinput >> gridy_input[cnt];
    cnt++;
  }
  
  if(z_coord){
    // Read the data from the file
    cnt = 0;
    for (k=0;k<nz;k++){
      std::getline(iss, line);
      std::istringstream sinput(line);
      sinput >> gridz_input[cnt];
      // gridz_input[cnt] = gridz_input[cnt]; //convert to [cm]
      cnt++;
    }
  }

  // Read the data from the file
  cnt = 0;
  //while (std::getline(iss, line)) {
  for (n=0;n<nscal;n++){
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
        for (i=0;i<nx;i++){
          std::getline(iss, line);
          std::istringstream sinput(line);
          //read only the specified number of points in x defined by the array ix_dir_limit
          if(i >= ix_dir_limit[0] and i < ix_dir_limit[1]){
            sinput >> data[cnt]; 
            cnt++;
          }
        }
      }
    }
  }
  //if(cnt != PeleC::h_prob_parm_device->nscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny){
  //  printf("Number of lines = %i. Specified in the header = %i \n", cnt,PeleC::h_prob_parm_device->nscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny);
  //  amrex::Abort("Number of lines in file differ from specified in header");
  //}
};