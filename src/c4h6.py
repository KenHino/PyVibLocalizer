from pyvib import simulator, read_minfo, write_minfo #geometry, freq, displacement vector


#geom, freq, disp = read_minfo("minfo_sample/butadiene.minfo", use_trans=True, use_rot=True)
geom, freq, disp = read_minfo("minfo_sample/butadiene.minfo")




#vertical_mode = [0,3,4,6,7,9,10]
#sim_vertical = simulator(geom, [freq[m] for m in vertical_mode], [disp[m] for m in  vertical_mode])
#sim_vertical.localize(option='Boys')


#horizontal_mode = [1,2,5,8,11,12,13,14,15,16,17]
#horizontal_mode = [11,12,13,14,15,16,17]
#sim_horizontal = simulator(geom, [freq[m] for m in horizontal_mode], [disp[m] for m in  horizontal_mode])
sim_horizontal = simulator(geom, freq[11:], disp[11:])
sim_horizontal.localize(option='Boys', window=200)
sim_horizontal.visualize()


#sim_CH_strech = simulator(geom, freq[18:], disp[18:])
#sim_CH_strech.localize(option='Boys')

#sim.visualize()
#sim.localize(option='Pipek-Mezy')
#sim = simulator(geom, sim_vertical.freq + sim_horizontal.freq + sim_CH_strech.freq, sim_vertical.disp + sim_horizontal.disp + sim_CH_strech.disp)
#sim = simulator(geom, sim_horizontal.freq + sim_CH_strech.freq, sim_horizontal.disp + sim_CH_strech.disp)
#sim.visualize()
#write_minfo("butadiene.out", sim.freq, sim.disp)
