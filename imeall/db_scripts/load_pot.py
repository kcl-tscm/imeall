from quippy import AtomsReader, Potential
r_scale =1.08
eam_pot ='/home/lambert/pymodules/imeall/imeall/potentials/PotBH.xml'
pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
