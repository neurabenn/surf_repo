#!/usr/bin/env python

import os,sys,argparse,shutil

import sys
import xml.etree.ElementTree as ET
import subprocess as sp
from pathlib import Path
fsldir=os.getenv('FSLDIR')
FS_dir=os.getenv('FREESURFER_HOME')


### get arguments and parse them
parser = argparse.ArgumentParser(description='Map a Volumetric Atlas to a Surface in Freesurfer space.',usage='atlas_2_surf.py -r <T1 template aligned to atlas> -a <atlas image file>  -l <xml atlas file> -s <subject to atlas to surface>',epilog=("Example usage: "+"atlas_2_surf.py -r templ.nii.gz -x atlas.xml -s "),add_help=False)
if len(sys.argv) < 2:
	err="{cmd}"
	err=err.format(cmd=sys.argv[0])
	sp.call(["python",err,"--help"])
	sys.exit(1)

req_grp = parser.add_argument_group(title='Required arguments')
req_grp.add_argument('-r','--ref',type=str,metavar='',required=True,help='Template Image Associated with atlas labels')
req_grp.add_argument('-a','--atlas',type=str,metavar='',required=True,help='Atlas image file. Typically a set of binary labels')
req_grp.add_argument('-l','--labels',type=str,metavar='',required=True,help='Associated XML file to Atlas Image')
req_grp.add_argument('-s','--surf',type=str,metavar='',required=True,help='Statistical or stat image to sample for output PNG file')


#parser.add_argument('-a','--anat',type=str,metavar='',required=True,help='Anatomical Image aligned to mask space. Default is to sample images here for png output')
#parser.add_argument('-m','--mask',type=str,metavar='',required=True,help='Statistical or stat image to sample for output PNG file')
args=parser.parse_args()

#get those arguments into variables
temp_vol=args.ref
atl=args.atlas
lab=args.labels
surf=args.surf
cwd = os.getcwd()

#determine the midline for splitting the left and right hemisheres using the anatomic template the atlas is aligned too
mid_cmd="{fsl}/bin/fslstats {img} -C".format(fsl=fsldir,img=temp_vol)
mid_line=sp.run(mid_cmd.split(),stdout=sp.PIPE)
mid_line=mid_line.stdout.decode('utf-8')
mid_line=mid_line.split('.')
mid_line=int(mid_line[0])

print("The line for splitting the left and right hemispheres is {middle}".format(middle=mid_line))


#### conform the atlas into freesurfer space.
#### future generations can be nonlinearly warped to do indivdual labels. for now though this is from template to template only
FS_atl="FS_{ext}".format(ext=atl.split('/')[-1])

conf_cmd="{FS}/bin/mri_convert {img} {FS_atl} -rt nearest -nc --conform".format(FS=FS_dir,img=atl,FS_atl=FS_atl)
sp.run(conf_cmd.split())

####now we are ready to parse the xml and map the volumetric labels to their respective surfaces
## parse the xml atlas file


def map_mask_to_surface(img,idx,name,hemi,srf):
	#set subjects DIR to where you are
	os.environ['SUBJECTS_DIR']=cwd
	tmp="{dir}/tmp_label".format(dir=cwd)
	print(tmp)
	os.mkdir(tmp)
	print(os.getenv('SUBJECTS_DIR'))
	
	cmd1="{FS}/bin/mri_vol2label --i {img}  --id {idx} --l {tmp}/{hemi}.{name}.label".format(FS=FS_dir,img=img,idx=idx,tmp=tmp,hemi=hemi,name=name)
	cmd2="{FS}/bin/mri_label2vol --label {tmp}/{hemi}.{name}.label --temp {img} --o {tmp}/{hemi}.{name}_vol.mgh --identity".format(FS=FS_dir,img=img,idx=idx,tmp=tmp,hemi=hemi,name=name)
	cmd3="{FS}/bin/mri_vol2surf --mov {tmp}/{hemi}.{name}_vol.mgh --ref orig.mgz --hemi {hemi} --o {tmp}/{hemi}.{name}_srf.mgh  --trgsubject {srf} --regheader {srf}".format(FS=FS_dir,img=img,idx=idx,tmp=tmp,hemi=hemi,name=name,srf=srf)
	cmd4="{FS}/bin/mri_vol2label  --i {tmp}/{hemi}.{name}_srf.mgh   --id 1  --l avg_surf/label/{hemi}.{name}  --surf {srf} {hemi}".format(FS=FS_dir,img=img,idx=idx,tmp=tmp,hemi=hemi,name=name,srf=srf)
	sp.run(cmd1.split())
	sp.run(cmd2.split())
	sp.run(cmd3.split())
	sp.run(cmd4.split())


	shutil.rmtree(tmp)

lab_tree=ET.parse(lab)
lab_root=lab_tree.getroot()
for i in lab_root:
	for reg in i:
		data=reg.attrib
		if len(data) == 0:
			pass
		else:
			idx=data['index']
			x=int(data['x'])
			name=reg.text.replace(" ","_")
			if x < mid_line:
				hemi='rh'
				# print("the index of {name} is {idx}, with an x value of {x} it pertains to the {hemi}".format(name=name,idx=idx,x=x,hemi=hemi))
				### insert functon to register to surface here
				map_mask_to_surface(atl,idx,name,hemi,surf)
			else:
				hemi='lh'
				
				# print("the index of {name} is {idx}, with an x value of {x} it pertains to the {hemi}".format(name=name,idx=idx,x=x,hemi=hemi))
				### insert functon to register to surface here
				map_mask_to_surface(atl,idx,name,hemi,surf)



























