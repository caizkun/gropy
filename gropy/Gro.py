class Gro:
    """
        the central class in gropy
    """
    # -- constructor(s) --
    def __init__(self, 
        system_name=None, num_of_atoms=None, 
        residue_id=None, residue_name=None,
        atom_name=None, atom_id=None, 
        x=None, y=None, z=None, 
        v_x=None, v_y=None, v_z=None, 
        box=None):
        """
            wrap the contents in a GROMACS gro file in a class
        """
        self.system_name  = system_name or 'This is a Gro!'
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
    # not mandatory in python

    # -- file i/o --
    def read_gro_file(self, file_name):
        """
            read a gro file and store information in a Gro object
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


    def write_gro_file(self, file_name):
        """
            write a gro file based on a Gro object
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


    # -- preservative operations --
    def rename_atoms(self, old_atom_names, new_atom_names):
        """
            rename atoms from the old_atom_names to the new_atom_names
        """
        assert (len(old_atom_names) == len(new_atom_names)), "old_atom_names doesn't have the same length as new_atom_names"
        for i_atom in xrange(self.num_of_atoms):
            for i_name in xrange(len(old_atom_names)):
                if self.atom_name[i_atom] == old_atom_names[i_name]:
                    self.atom_name[i_atom] = new_atom_names[i_name]
                    break

    # TODO: may add flexibility to rename atoms with specific residue_names

    def rename_residues(self, old_residue_names, new_residue_names):
        """
            rename residues with an old_residue_name to a new_residue_name
        """
        assert (len(old_residue_names) == len(new_residue_names)), "old_residue_names doesn't have the same length as new_residue_names"
        for i_atom in xrange(self.num_of_atoms):
            for i_name in xrange(len(old_residue_names)):
                if self.residue_name[i_atom] == old_residue_names[i_name]:
                    self.residue_name[i_atom] = new_residue_names[i_name]
                    break


    def renumber_atoms(self):
        """
            renumber residue_id and atom_id starting from 1; the original composition of each resdiue is maintained
        """
        last_old_residue_id = -1                # use negative num to avoid coincidence
        last_old_residue_name = 'to_be_defined'
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


    def replace_atom_entry(self, i_atom, another_gro_object, j_atom):
        """
            replace the i-th atom of the current gro object with the j-th atom of another gro object
        """
        self.residue_id[i_atom]   = another_gro_object.residue_id[j_atom]
        self.residue_name[i_atom] = another_gro_object.residue_name[j_atom]
        self.atom_name[i_atom]    = another_gro_object.atom_name[j_atom]
        self.atom_id[i_atom]      = another_gro_object.atom_id[j_atom]
        self.x[i_atom]            = another_gro_object.x[j_atom]
        self.y[i_atom]            = another_gro_object.y[j_atom]
        self.z[i_atom]            = another_gro_object.z[j_atom]
        self.v_x[i_atom]          = another_gro_object.v_x[j_atom]
        self.v_y[i_atom]          = another_gro_object.v_y[j_atom]
        self.v_z[i_atom]          = another_gro_object.v_z[j_atom]


    def sort_residues(self, residue_name_list):
        """
           sort residues in the provided order, attaching other unspecified residues to the end
        """
        sorted_gro = Gro()
        # copy provided to another
        for residue_name in residue_name_list:
            sorted_gro.copy_residues(self, [residue_name])
        # remove provided
        self.remove_residues(residue_name_list)
        # copy unprovided to another
        for i_atom in xrange(self.num_of_atoms):
            sorted_gro.copy_atom_entry(self, i_atom)
        # delete unprovided
        for i_atom in xrange(self.num_of_atoms):
            self.remove_atom_entry(0)
        # copy all back from another
        for i_atom in xrange(sorted_gro.num_of_atoms):
            self.copy_atom_entry(sorted_gro, i_atom)


    # -- additive operations --
    def copy_atom_entry(self, another_gro_object, i_atom):
        """
            copy the i-th atom entry from another gro object and append to the end of current gro object
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
    
    
    def copy_residue_entry(self, another_gro_object, residue_id, residue_name):
        """
            copy atoms of the specified residue from another gro object and append to the end of current gro object
        """
        for i_atom in xrange(another_gro_object.num_of_atoms):
            if another_gro_object.residue_id[i_atom] == residue_id and another_gro_object.residue_name[i_atom] == residue_name:
                self.copy_atom_entry(another_gro_object, i_atom)


    def copy_atoms(self, another_gro_object, atom_name_list):
        """
            copy atoms with the provided atom names from another gro object and append to the end of current gro object
        """
        for i_atom in xrange(another_gro_object.num_of_atoms):
            for atom_name in atom_name_list:
                if another_gro_object.atom_name[i_atom] == atom_name:
                    self.copy_atom_entry(another_gro_object, i_atom)
                    break


    def copy_residues(self, another_gro_object, residue_name_list):
        """
            copy atoms with the provided residue names from another gro object and append to the end of current gro object
        """
        for i_atom in xrange(another_gro_object.num_of_atoms):
            for residue_name in residue_name_list:
                if another_gro_object.residue_name[i_atom] == residue_name:
                    self.copy_atom_entry(another_gro_object, i_atom)
                    break

    # TODO: may add to copy atoms with the providied atom names and residue names

    # Python doesn't allow the overloading of assignment operator "=";
    # In python, copying an object is often achived by utilizing the copy and deepcopy functions in the copy module.


    # -- subtractive operations --
    def remove_atom_entry(self, i_atom):
        """
            remove the i-th atom entry from current gro object
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


    def remove_residue_entry(self, residue_id, residue_name):
        """
            remove atoms of the specified residue
        """
        atom_indice_to_be_removed = []
        for i_atom in xrange(self.num_of_atoms):
            if self.residue_id[i_atom] == residue_id and self.residue_name[i_atom] == residue_name:
                atom_indice_to_be_removed.append(i_atom)                        # save indice first; direct removal would shrink the atom list

        num_of_atoms_to_be_removed = len(atom_indice_to_be_removed)
        for i_atom in xrange(num_of_atoms_to_be_removed):
            self.remove_atom_entry(atom_indice_to_be_removed[i_atom] - i_atom)  # shift atom indice to match the shrinkage of atom list


    def remove_atoms(self, atom_name_list):
        """
            remove atoms with the provided atom names
        """
        atom_indice_to_be_removed = []
        for i_atom in xrange(self.num_of_atoms):
            for atom_name in atom_name_list:
                if self.atom_name[i_atom] == atom_name:
                    atom_indice_to_be_removed.append(i_atom)
                    break
        num_of_atoms_to_be_removed = len(atom_indice_to_be_removed)
        for i_atom in xrange(num_of_atoms_to_be_removed):
            self.remove_atom_entry(atom_indice_to_be_removed[i_atom] - i_atom)  # shift atom indice to match the shrinkage of atom list


    def remove_residues(self, residue_name_list):
        """
            remove atoms with the provided residue names
        """
        atom_indice_to_be_removed = []
        for i_atom in xrange(self.num_of_atoms):
            for residue_name in residue_name_list:
                if self.residue_name[i_atom] == residue_name:
                    atom_indice_to_be_removed.append(i_atom)
                    break
        num_of_atoms_to_be_removed = len(atom_indice_to_be_removed)
        for i_atom in xrange(num_of_atoms_to_be_removed):
            self.remove_atom_entry(atom_indice_to_be_removed[i_atom] - i_atom)  # shift atom indice to match the shrinkage of atom list

    # TODO: may add to copy atoms with the providied atom names and residue names
