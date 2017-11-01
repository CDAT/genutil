from __future__ import print_function
import cdms2,genutil,os,sys
import unittest
import cdat_info
from genutil.statistics import StatisticsError
import numpy

class GENUTIL(unittest.TestCase):
    def setUp(self):
        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'clt.nc'))
        self.clt=f('clt')
        f.close()
        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'tas_cru_1979.nc'))
        self.tas=f('tas',time=slice(0,2))
        f.close()
        self.varResult = [619.3521759227226, 614.4491584168575, 588.1676093151207, 603.488815605165,
             648.5501168851572, 686.6054391292984, 701.5697076183785, 683.1722831771162,
              586.0022620934736, 582.1263027572609, 588.1230361494828, 617.844886887937,
               630.6878165232655, 630.0147623402455, 599.205625962692, 612.3671105076764,
                668.1353673715317, 709.9454658867054, 685.757266660344, 653.6661984981,
                 597.6718046928563, 538.2119547320317, 566.0159036732502, 612.8417453219153,
                  604.3785792308794, 617.5273929283923, 587.2523217329198, 612.2952451498608,
                   683.6324840170314, 721.7339270677919, 698.468637784132, 661.8158665533944,
                    575.3558802853665, 517.3843882873233, 579.4092887307754, 598.7339774672442,
                     617.011050713242, 618.9126321246761, 607.8986935607807, 619.6528777012809,
                      660.7785876679729, 707.3360504818712, 712.3877705490231, 646.2624139240238,
                       571.8818826451019, 566.9679871463228, 585.1198459882157, 623.8608207815935,
                        606.0901297851185, 632.387444382224, 621.7858380580385, 544.1908976793201,
                         613.508170281664, 668.5471178042499, 685.8647183767764, 683.8944732305736,
                          623.2817871937893, 572.0952313656045, 580.0422629685613, 632.0045397153427,
                           642.1948460458397, 642.4829257450631, 611.7783510272177, 606.072801202646,
                            656.1744779521672, 674.1221654887953, 704.1429322473252, 672.6335331178888,
                             623.2337248842003, 572.2441559781086, 612.5786757628019, 613.9894036696215,
                              633.7451517234747, 633.8605771060166, 633.2874412467597, 590.1857030748076,
                               665.619360595367, 725.3431232665431, 706.2409202309975, 679.0585511805239,
                                586.2206116289849, 544.7416975764929, 591.2834795462458, 639.7088714354439,
                                 651.2303685776834, 629.179528634222, 608.2348369843769, 604.5912990703287,
                                  672.3113243166489, 688.833001363026, 706.0611142859545, 685.9627800771833,
                                   620.9449931279593, 569.3150375391932, 616.4511545370198, 613.7642290855487,
                                    641.1817429213651, 599.2781484489318, 581.0695106700597, 627.3442135194779,
                                     665.9524887789161, 687.0931824802278, 695.2664060325218, 704.650835729694,
                                      575.3011978645626, 550.4301110373642, 608.6312986716732, 652.8933354700406,
                                       634.8464859898008, 600.5446776647464, 589.5394153523209, 582.2753179660784,
                                        627.6516102875161, 692.4160407146479, 700.4200948975626, 656.1743898418398,
                                         591.9393094998022, 582.1433460781363, 605.8020207962994, 609.6489688600152]


    def testChecker_1(self):
        axistest=[0,[1,2],'xy','12','(longitude)(latitude)']
        for ax in axistest:
            print('testing:',ax)
            v=genutil.statistics.variance(self.tas,axis=ax)

        ax='(dsffd)(latitude)'
        with self.assertRaises(StatisticsError):
            v=genutil.statistics.variance(self.tas,axis=ax)

    def testChecker_2(self):
        print("####################### Test 2 ################################")
        v=genutil.statistics.variance(self.tas,axis='012',weights=['unweighted','weighted','unweighted'])
        self.assertEqual(v.ndim,0)
        self.assertTrue(numpy.allclose(float(v),218.880632191))
        v=genutil.statistics.variance(self.tas,axis='12',weights='weighted')
        self.assertEqual(v.ndim,1)
        self.assertTrue(numpy.allclose(v,[220.10595800589059, 217.5567693498757]))

    def testChecker_3(self):
        print("####################### Test 3 ################################")
        v=genutil.statistics.covariance(self.clt,self.clt,axis='(latitude)(longitude)',weights=['generate','generate'])
        self.assertEqual(v.ndim,1)
        self.assertTrue(numpy.allclose(v,self.varResult))

    def testChecker_4(self):
        print("####################### Test 4 ################################")
        a=genutil.statistics.variance(self.clt, axis='tyx',weights=['equal', 'generate', 'equal'])
        self.assertEqual(a.ndim,0)
        self.assertTrue(numpy.allclose(a,630.549411264))
        a=genutil.statistics.variance(self.clt, axis='yx',weights=['generate', 'equal'])
        self.assertEqual(a.ndim,1)
        self.assertTrue(numpy.allclose(a,self.varResult))
        a=genutil.statistics.variance(self.clt, axis='tx',weights=['equal', 'equal'])
        self.assertTrue(numpy.allclose(a,[107.68428802490234, 492.9377746582031, 425.12725830078125, 406.87811279296875,
             294.3303527832031, 328.2757568359375, 457.69036865234375, 257.73638916015625,
              78.1517105102539, 12.577239036560059, 25.770889282226562, 53.896240234375,
               107.05201721191406, 160.54620361328125, 282.6634826660156, 391.7831115722656,
                436.4649353027344, 441.2729797363281, 421.15972900390625, 373.83819580078125,
                 369.6070556640625, 309.7768249511719, 327.92578125, 298.1361083984375,
                  309.5687561035156, 395.45379638671875, 489.1477355957031, 627.8700561523438,
                   748.49560546875, 788.53662109375, 830.0309448242188, 743.3065795898438,
                    684.14013671875, 747.7320556640625, 776.3470458984375, 755.2000122070312,
                     631.7830810546875, 522.89404296875, 343.2276306152344, 257.6018371582031,
                      305.44769287109375, 379.7539367675781, 404.60992431640625, 420.9342346191406,
                       497.81341552734375, 688.7490234375]))
        a=genutil.statistics.variance(self.clt, axis='t',weights=['equal'])
        self.assertEqual(a.shape,(46,72))


    def testChecker_5(self):
        print("####################### Test 5 ################################")
        a=genutil.statistics.laggedcovariance(self.clt,self.clt,lag=4,axis='(time)(longitude)')
        self.assertEqual(a.shape,(5,46))

