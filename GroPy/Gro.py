class GroSystem:
    """
    wrap the contents in a GROMACS gro file as a class
    """
    # -- constructor(s) --
    def __init__(self, 
        system_name=None, num_of_atoms=None, 
        residue_id=None, residue_name=None,
        atom_name=None, atom_id=None, 
        x=None, y=None, z=None, 
        v_x=None, v_y=None, v_z=None, 
        box=None):
        self.system_name  = system_name or 'This is a Gro System'
        self.num_of_atoms = num_of_atoms or 0
        self.residue_id   = residue_id or []
        self.residue_name = residue_name or []
        self.atom_name    = atom_name or []
        self.atom_id      = atom_id or []
        self.x            = x or []
        self.y            = y or []
        self.z            = z or []
        self.v_x          = v_x or []
        self.v_y          = v_y or []
        self.v_z          = v_z or []
        self.box          = box or [0.0, 0.0, 0.0]


    # -- deconstructor --
    # not necessary in python
    

    # -- file i/o --
    def read_gro_file(self, file_name):
        """
        read a gro file and store information in a GroSystem object
        """
        with open(file_name, 'r') as file_id:
            for i_line, line in enumerate(file_id):
                line = line.replace('\n', '')

                if i_line == 0:
                    self.system_name = line
                    continue

                if i_line == 1:
                    self.num_of_atoms = int(line)
                    final_line_of_atoms = self.num_of_atoms + 1
                    continue

                if i_line <= final_line_of_atoms:
                    # store atom information
                    self.residue_id.append(int(line[0:5]))
                    self.residue_name.append(line[5:10].strip())    # remove leading spaces
                    self.atom_name.append(line[10:15].strip())      # remove leading spaces
                    self.atom_id.append(int(line[15:20]))
                    self.x.append(float(line[20:28]))
                    self.y.append(float(line[28:36]))
                    self.z.append(float(line[36:44]))
                    if len(line) > 44:
                        self.v_x.append(float(line[44:52]))
                        self.v_y.append(float(line[52:60]))
                        self.v_z.append(float(line[60:68]))
                    else:
                        self.v_x.append(0.0)
                        self.v_y.append(0.0)
                        self.v_z.append(0.0)
                else:
                    self.box = line.split()
                    self.box = [float(box_size) for box_size in self.box]

        file_id.close()


    def write_gro_file(self, file_name):
        """
        write a gro file based on a GroSystem object
        """
        with open(file_name, 'w') as file_id:
            file_id.write("%s\n" % self.system_name)
            file_id.write(" %d\n" % self.num_of_atoms)
            
            if self.v_x != []:
                for i in xrange(self.num_of_atoms):
                    file_id.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (self.residue_id[i], self.residue_name[i], 
                        self.atom_name[i], self.atom_id[i], self.x[i], self.y[i], self.z[i], self.v_x[i], self.v_y[i], self.v_z[i]))
            else:
                for i in xrange(self.num_of_atoms):
                    file_id.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (self.residue_id[i], self.residue_name[i], 
                        self.atom_name[i], self.atom_id[i], self.x[i], self.y[i], self.z[i], 0.0, 0.0, 0.0))
            
            file_id.write("%10.5f%10.5f%10.5f\n" % (self.box[0], self.box[1], self.box[2]))

        file_id.close()


    # -- preservative operations --
    def rename_atoms(self, old_atom_name, new_atom_name):
        """
        rename atoms with an old_atom_name to a new_atom_name
        """
        for i_atom in xrange(self.num_of_atoms):
            if self.atom_name[i_atom] == old_atom_name:
                self.atom_name[i_atom] = new_atom_name
                    
    
    def rename_residues(self, old_residue_name, new_residue_name):
        """
        rename residues with an old_residue_name to a new_residue_name
        """
        for i_atom in xrange(self.num_of_atoms):
            if self.residue_name[i_atom] == old_residue_name:
                self.residue_name[i_atom] = new_residue_name


    def renumber_atoms(self):
        """
        renumber residue_id and atom_id starting from 1; original constitution of each resdiue is maintained
        """
        last_old_residue_id = -1                # use negative num to avoid coincidence
        last_old_residue_name = 'ToBeDefined'
        last_new_resiue_id = 0
        for i_atom in xrange(self.num_of_atoms):
            self.atom_id[i_atom] = i_atom + 1   # starting from 1
            if self.residue_id[i_atom] == last_residue_id and self.residue_name[i_atom] == last_residue_name:
                self.residue_id[i_atom] = last_new_resiue_id
            else:
                last_old_residue_id = self.residue_id[i_atom]
                last_old_residue_name = self.residue_name[i_atom]
                self.residue_id[i_atom] = last_new_resiue_id + 1
                last_new_resiue_id += 1

## TODO: finish sort

#    def sort_by_atom_name(self, atom_name_list):
#        """
#        sort atoms by their atom names in the order of provided atom_name_list, moving unspecified atoms to the end
#        """
#        system_of_sorted_atoms = GroSystem()
#        for specified_atom_name in atom_name_list:
#            for i_atom in xrange(self.num_of_atoms):
#                if self.atom_name[i_atom] == specified_atom_name:
#                    system_of_sorted_atoms.copy_atom_entry(self, i_atom)


#    def sort_by_residue_name(self, residue_name_list):
#        """
#        sort atoms by their residue names in the order of provided residue_name_list, , moving unspecified atoms to the end
#        """


    # -- additive operations --
    def copy_atom_entry(self, another_gro_object, i_atom):
        """
        copy one atom entry from another gro object and append to the end of current gro object
        """
        self.num_of_atoms += 1
        self.residue_id.append(another_gro_object.residue_id[i_atom]);
        self.residue_name.append(another_gro_object.residue_name[i_atom]);
        self.atom_name.append(another_gro_object.atom_name[i_atom]);
        self.atom_id.append(another_gro_object.atom_id[i_atom]);
        self.x.append(another_gro_object.x[i_atom]);
        self.y.append(another_gro_object.y[i_atom]);
        self.z.append(another_gro_object.z[i_atom]);
        self.v_x.append(another_gro_object.v_x[i_atom]);
        self.v_y.append(another_gro_object.v_y[i_atom]);
        self.v_z.append(another_gro_object.v_z[i_atom]);
    
    
    def copy_residue_entries(self, another_gro_object, residue_id, residue_name):
        """
        copy atoms within specified residue from another gro object and append to the end of current gro object
        """
        for i_atom in xrange(another_gro_object.num_of_atoms):
            if another_gro_object.residue_id[i_atom] == residue_id and another_gro_object.residue_name[i_atom] == residue_name:
                self.copy_atom_entry(another_gro_object, i_atom)


    def copy_atoms_by_name(self, another_gro_object, atom_name):
        """
        copy atoms with the provided atom name from another gro object and append to the end of current gro object
        """
        for i_atom in xrange(another_gro_object.num_of_atoms):
            if another_gro_object.atom_name[i_atom] == atom_name:
                self.copy_atom_entry(another_gro_object, i_atom)


# TODO: add merge systems; add copy_by_atom_name


    # -- subtractive operations --
    def remove_atom_entry(self, i_atom):
        """
        remove i-th atom from current gro object
        """
        self.num_of_atoms -= 1
        del self.residue_id[i_atom]
        del self.residue_name[i_atom]
        del self.atom_name[i_atom]
        del self.atom_id[i_atom]
        del self.x[i_atom]
        del self.y[i_atom]
        del self.z[i_atom]
        del self.v_x[i_atom]
        del self.v_y[i_atom]
        del self.v_z[i_atom]


    def remove_residue_entries(self, residue_id, residue_name):
        """
        remove atoms within specified residue from current gro object
        """
        atom_indice_to_be_removed = []
        for i_atom in xrange(self.num_of_atoms):
            if self.residue_id[i_atom] == residue_id and self.residue_name[i_atom] == residue_name:
                atom_indice_to_be_removed.append(i_atom)                        # save indice first; direct removal would shrink the atom list

        num_of_atoms_to_be_removed = len(atom_indice_to_be_removed)
        for i_atom in xrange(num_of_atoms_to_be_removed):
            self.remove_atom_entry(atom_indice_to_be_removed[i_atom] - i_atom)  # shift atom indice to match the shrinkage of atom list


    def remove_atoms_by_name(self, atom_name):
        """
        remove atoms with the provided atom name
        """
        atom_indice_to_be_removed = []
        for i_atom in xrange(self.num_of_atoms):
            if self.atom_name[i_atom] == atom_name:
                atom_indice_to_be_removed.append(i_atom)

        num_of_atoms_to_be_removed = len(atom_indice_to_be_removed)
        for i_atom in xrange(num_of_atoms_to_be_removed):
            self.remove_atom_entry(atom_indice_to_be_removed[i_atom] - i_atom)  # shift atom indice to match the shrinkage of atom list


