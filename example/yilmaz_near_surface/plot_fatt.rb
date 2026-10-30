
iter = 100

opts = "-n1=21 -label1='Depth (m)' -label2='Horizontal Position (m)' -legend=1 -lloc=bottom -unit='P-wave Velocity (m/s)' -cmin=350 -cmax=5500 -size1=2.5 -size2=6 -lmtick=9 -ld=1000 -tick1d=5 -mtick1=4 -tick2d=10 -mtick2=9 -interp=Gaussian "

system "mkdir -p result"

#system "x_showmatrix -in=model/vp_init.bin #{opts} -out=./result/vp_0.png -title='Initial' &"

#system "cat pick/shot_1_*.bin pick/shot_10_* pick/shot_20_* pick/shot_30_* pick/shot_40_* >obs.bin "

#for i in 1..100
	
	#system "x_showmatrix -in=test_fatt/iteration_#{i}/model/updated_vp.bin #{opts} -title='Iteration #{i}' -out=./result/vp_#{i}.png "
	
	#dir = "test_fatt/iteration_#{i}/synthetic"
	#system "cat #{dir}/shot_1_*.bin #{dir}/shot_10_* #{dir}/shot_20_* #{dir}/shot_30_* #{dir}/shot_40_* >syn_#{i}.bin "
	#system "x_showgraph -in=obs.bin,syn_#{i}.bin -n1=240,240 -x2beg=0 -x2end=0.05 -tick2d=0.01 -mtick2=9 -label2='Time (s)' -label1='Trace Number' -tick1d=25 -mtick1=4 -linewidth=2,2 -linecolor=b,r -reverse2=1 -size1=8 -size2=2.5 -plotlabelloc=lower_right -plotlabel='Observed':'Iteration #{i}' -out=./result/data_#{i}.png & "
		
#end

system "ffmpeg -framerate 10 -i ./result/vp_%d.png -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -c:v libx264 -pix_fmt yuv420p ./vp.mp4"
system "ffmpeg -framerate 10 -i ./result/data_%d.png -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -c:v libx264 -pix_fmt yuv420p ./data.mp4"
