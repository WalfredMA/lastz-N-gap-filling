#!/usr/bin/python

import os
import pandas as pd
import re
import collections as cl


opts,args=getopt.getopt(sys.argv[1:],"i:o:r:")
for op, value in opts:
	if op=='-i':
		querypath=value
	if op=='-o':
		outputpath=value
	if op=='-r':
		refpath=value

size0=100
alignment0=2000
insert0=50
match=98
match0=49	


def distractgap(x,nongapcount=9):
	
	t=x.replace('-','N')
	a1=re.finditer(r'N+',t)
	a1=[a.span() for a in a1]
	if len(a1) == 1:
		return [x[a1[0][0]:a1[0][1]]],[a1[0][0]],[a1[0][1]]
	if len(a1) == 0:
		return [],[],[]
	a2=[[a1[i][1]-a1[i][0],a1[i+1][0]-a1[i][1],a1[i+1][1]-a1[i+1][0]] for i in range(0,len(a1)-1)]
	nonignore=[]
	for i in range(len(a2)/2):
		m1=i
		m2=len(a2)-1-i
		s0=(a2[m1][1])*nongapcount
		s1=(a2[m2][1])*nongapcount
		if s0 > a2[m1][0] or s0 > a2[m1][2]:
			nonignore.append(m1)
			if s1 > a2[m2][0] or s1 > a2[m2][2]:
				nonignore.append(m2)
				continue
			a2[m2-1][2]+=a2[m2][2]-s0
			continue
		if s1 > a2[m2][0] or s1 > a2[m2][2]:
			nonignore.append(m2)
			a2[m2-1][2]+=a2[m2][2]-s0
			continue
		a2[m2-1][2]+=a2[m2][2]-s0
		a2[m1+1][0]+=a2[m1][0]-s1
	if len(a2)/2*2 != len(a2):
		m1=len(a2)-len(a2)/2-1
		s0=(a2[m1][1])*nongapcount
		if s0 > a2[m1][0] or s0 > a2[m1][2]:
			nonignore.append(m1)
	nonignore.sort()
	p1=[a1[a][1] for a in nonignore]
	p1.append(a1[-1][1])
	p0=[a1[a+1][0] for a in nonignore]
	p0.insert(0,a1[0][0])
	
	
	return p0, p1


def blastncheck(query,ref,lr):
	
	def align_combine(qstart,qend,rchr,rstd,rstart,rend):
		
		groups=cl.defaultdict(list)
		
		for i,chr0 in enumerate(rchr):
			
			groups[chr0+':'+rstd[i]].append(i)
					
		alignments=[]
		for key in groups.keys():
			
			g_index=groups[key]
			
			g_rstart=[rstart[i] for i in g_index]
			
			g_rend=[rend[i] for i in g_index]
			
			alignments_index0=eachgroup_combine(g_index,g_rstart,g_rend,10000)
			
			
			for each_index in alignments_index0:

				
				each_qstart=[qstart[i] for i in each_index]
				
				each_qend=[qend[i] for i in each_index]
				
				
				each_index_q=eachgroup_combine(range(len(each_qstart)),each_qstart,each_qend,20)
				
				each_qstart=[min([each_qstart[x] for x in  index0]) for index0 in each_index_q]
				
				each_qend=[max([each_qend[x] for x in  index0]) for index0 in each_index_q]
				
				chr0=key.split(':')[0]
				
				std0=key.split(':')[1]
				
				if std0=='+':
							
					each_rstart=min([rstart[i] for i in each_index])
					
					each_rend=max([rend[i] for i in each_index])
					
				
				else:
					
					each_rstart=max([rstart[i] for i in each_index])
					
					each_rend=min([rend[i] for i in each_index])
					
				size0=sum([abs(each_qend[i]-each_qstart[i]) for i in xrange(len(each_qend))])
				
				alignments.append([chr0,std0,size0,each_qstart,each_qend,each_rstart,each_rend])
		
		return alignments

	
	
	os.system('blastn -query {:s} -db {:s} -outfmt 10 -out {:s} -num_threads 1 -max_target_seqs 150'.format(query, ref, query+'out'))	

	try:
		t=pd.read_csv(query+'out',sep=',',header=None)
	except:
		return  [],0

	if len(t)<1:

		return [],0


	qstart,qend,rchr,rstart,rend=list(t[6]),list(t[7]),list(t[1]),list(t[8]),list(t[9])

	rstd=['+' if rstart[i]<rend[i] else '-' for i in xrange(len(t))]

	alignments=align_combine(qstart,qend,rchr,rstd,rstart,rend)

	t=pd.DataFrame.from_records(alignments)

	nmatch=list(t[2])

	index=[i for i in xrange(len(t)) if nmatch[i]>(lr)*0.8]

	t=t.iloc[index]


	if len(t)<1:

		return []	


	results=[]
	for i in xrange(len(t)):


		contig=list(t[0])[i]

		size0=lr

		std=list(t[1])[i]
		
		qstart=min(list(t[3])[i])
		
		qend=max(list(t[4])[i])
		
		score=100*float(list(t[2])[i])/(lr)

		strand=list(t[1])[i]
		
		rstart1=int(list(t[5])[i])
		
		rend1=int(list(t[6])[i])

				
		results.append([contig,strand,rstart1,rend1,qstart,size0-qend,score])
	
	
		
	
	return sorted(results, key=lambda x: x[-1])




def findgap(refpath,tempfolder,file0,twofiles,sizes):
	
	
	t_front,t_back=[],[]
	
	if twofiles[0]==1:
		
		frontfind=blastncheck(twofiles[0],refpath,sizes[0])
		
			
	if twofiles[1]==1:
		
		endfind=blastncheck(twofiles[1],refpath,sizes[1])
		
	allsamechre=[]
	for i,frontfind0 in enumerate(frontfind):

		chr0,std,s0,e0=frontfind0[0],frontfind0[1],frontfind0[2],frontfind0[3]
		
		samechre=[a for a in endfind if a[0]==chr0]

		dis=[]
		for samechre0 in samechre:

			cutstart=min(s0,e0,samechre0[2],samechre0[3])

			cutend=max(s0,e0, samechre0[2],samechre0[3])
				
			dis0=cutend-cutstart-abs(s0-e0)-abs(samechre0[2]-samechre0[3])

			dis.append(dis0)
				
		samechreindex=sorted([x for x in xrange(len(samechre)) if dis[x]<100000],key=lambda x : dis[x])

		samechre=[samechre[x] for x in samechreindex]	

		if len(samechre)==0:

			continue

		bothscore=frontfind0[6]+samechre[0][6]

		allsamechre.append([frontfind0,samechre,bothscore])
		
	if len(allsamechre)>0:
		
		allsamechre=sorted(allsamechre,key=lambda x: x[-1])
		
		allsamechre[0][0]=file0+allsamechre[0][0]
		
		allsamechre[1][0]=file0+allsamechre[1][0]
		
		pd.DataFrame.from_records(allsamechre[0][:2]).to_csv(tempfolder+'/allblast.csv', header=None,sep=',',index=False, mode='a')
	
	
	return allsamechre

	
def prepare(querypath,outputpath):
	
	def cutfasta(file0,tempfolder):
		
		filename=file0.split('/')[-1]
		
		with open(file0,mode='r') as f:
			
			reads=f.read().split('>')[1:]
		
		f.close()
		
		for i,read in enumerate(reads):
			
			outputfile0='{:s}/{:s}~{:d}'.format(tempfolder,filename,i)
			
			with open(outputfile0, mode='w') as f:
				f.write('>'+read)
			f.close()
	
	tempfolder=outputpath+'/temp00'
	
	try:
		os.mkdir(outputpath+'/')
	except:
		pass	
	
	
	try:
		os.mkdir(tempfolder+'/')
	except:
		print "warnning, temp00 exits"
		
	if os.path.isfile(querypath)
	
		allfiles=[querypath]
		
	else:
		
		allfiles=os.listdir(querypath)
	
	for file0 in allfiles:
		
		cutfasta(file0,tempfolder)
	
	return tempfolder

def organizefiles(tempfolder,outputpath):
	
	allfiles=os.listdir(tempfolder)
	
	original_infor=cl.defaultdict(list)
	
	for file0 in allfiles:
		
		prefix=file0.split('+')[0]
		
		original_infor[prefix].append(file0)
		
	
	for eachfile in original_infor.keys():
		
		allcontigs=sorted(original_infor[eachfile], key=lambda x: x.split('+')[-1])
		
		for contigfile in allcontigs:
			
			os.system('cat {:s} >> {:s}'.format(tempfolder+'/'+contigfile, outputpath+'/'+eachfile))


def correctloop(tempfolder,refpath):
		
	allfiles=os.listdir(tempfolder)
	
	for i,x in enumerate(allfiles):
		
		print('start NO.'+str(i)+'th turn\n\nstart calling lastz\n')
		
		print('start correction\n')
		
		file0=tempfolder+'/'+x
		
		original_name=x.split('~')[0]
		sort=x.split('~')[1]
		
		f=open(file0, mode='r')
		read=f.read()
		f.close()
		
		query=read.split('\n')
		title=query[0].strip()
		query=''.join(query[1:]).strip().replace('\n','')
		
		allgaps=zip(*distractgap(query))
		
		allcuts=[((max(0, gap[0]-1000), gap[0]),(gap[1], min(len(query), gap[1]+1000)  )) for gap in allgaps if abs(gap[1]-gap[0])>100]
		
		
		for i,cuts in enumerate(allcuts):
			
			twofiles=['','']
			
			sizes=[abs(cuts[0][0]-cuts[0][1]), abs(cuts[1][0]-cuts[1][1])]
			
			frontfile=file0+'~'+str(i)+'f'
			backfile=file0+'~'+str(i)+'b'
			
			if sizes[0]>100:
				
				twofiles[0]=frontfile
				with open(frontfile,mode='w') as f:
					f.write('>front\n'+query[cuts[0][0]:cuts[0][1]])
				f.close()
				
			if sizes[1]>100:
				
				twofiles[1]=backfile
				with open(backfile,mode='w') as f:
					f.write('>back\n'+query[cuts[1][0]:cuts[1][1]])
				f.close()
			
			if sum(twofiles)>0:
				
				results=findgap(refpath,tempfolder,file0,twofiles,sizes)
				
				for result in results:
		
					front_break_ref=result[0][3]
					front_break_query=cuts[0][1]-result[0][5]
				
					back_break_ref=result[1][2]
					back_break_query=cuts[1][0]+result[1][4]
				
					out=[original_name,sort,i,cuts[0][1],cuts[1][0],front_break_query,back_break_query,front_break_ref,back_break_ref]
				
					with open('allngaps_find.csv', mode='a') as f:
					
						f.write(','.join(out)+'\n')
					
					f.close()
		

		
	return 0



def main(querypath,outputpath,refpath):

	tempfolder=prepare(querypath,outputpath)

	correctloop(tempfolder,refpath)

	organizefiles(tempfolder,outputpath)

	try:
		os.system('rm -rf temp00')
	except:
		print('warning: unable to delect temporary folder temp00, please do manually\n')

	print 'finished program\n'
	

#main fuction to run step by step 

if __name__=='__name__':
	
	main(querypath,outputpath,refpath)