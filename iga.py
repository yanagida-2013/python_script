#!/usr/bin/env python3
# coding: utf-8

'''
井頭研究室データ解析用のPythonライブラリ
1. リビン用のクラス※
2. anampaを通したリストデータのゲート用クラス（めっちゃ遅い）
3. その他いろんな便利な関数

使うためにはいくつかのモジュールのインストールが必要です。
- numpy
- scipy
- matplotlib

windowsならここからダウンロード
- http://www.lfd.uci.edu/~gohlke/pythonlibs/

※ リビンに関して
ペレトロンのNaI(Tl)を使っている場合は問題無いですが、
特殊な検出器を使って校正式が一次間数で無い場合、
このプログラムは使えません。
'''

__author__ = "s.yanagida"
__version__ = "0.1.3"
__date__ = "29 Jan 2015"


# 外部ライブラリ
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
# デフォルトのライブラリ
import struct
import os
import sys
import subprocess


class Rebin:
    """
    A class for rebinning
    """

    def __init__(self, arg):
        '''
        arg : filename or numpy.ndarray, 2 or 3 column data
        '''
        if type(arg) == str:
            try:
                self.data = np.loadtxt(arg)
                col, row = self.data.shape
                if row == 2:
                    self.__err = False
                elif row == 3:
                    self.__err = True
                else:
                    print("# of rows must be 2 or 3.")
                    sys.exit()
            except:
                print("File load Error")
                sys.exit()
        elif type(arg) == np.ndarray:
            col, row = arg.shape
            if row == 2:
                self.data = arg
                self.__err = False
            elif row == 3:
                self.data = arg
                self.__err = True
            else:
                print("# of rows must be 2 or 3.")
                sys.exit()
        else:
            print("Input error")
            sys.exit()

    def all(self):
        '''
        retrun all-data as numpy.ndarray
        '''
        return self.data

    def x(self):
        '''
        return x-axis
        '''
        return self.data[:, 0]

    def y(self):
        '''
        return y-axis
        '''
        return self.data[:, 1]

    def z(self):
        '''
        return z-axis
        '''
        if self.__err:
            return self.data[:, 2]
        else:
            return []

    def show(self):
        """
        Shows graph with matplotlib.
        Error bars will apper if z-axis exists.
        """
        if self.__err:
            plt.errorbar(self.x(), self.y(), yerr=self.z(),  marker='.')
        else:
            plt.scatter(self.x(), self.y(), marker='.')
        plt.show()

    def save(self, filename):
        """
        save as filename
        """
        np.savetxt(filename, self.data, fmt="%6.5e", delimiter="\t")

    def rebin(self, min, max, step):
        """
        Rebining function.
        min : minimum x-value
        max : miximum x-value
        step : step of x-value
        """

        # 古い値を取得
        old_x = self.x()
        old_y = self.y()
        old_z = self.z()
        # 新しい値を初期化
        if step < 1E-5:
            print("---------------------------------------------")
            print("Caution! Floating-point arithmetic may occur.")
            print("---------------------------------------------")
        new_x = np.arange(min, max + step*1E-5, step)
        mid_x = np.arange(min-step/2, max+step/2+step*1E-5, step)  # 切れ目
        new_y = np.zeros(len(new_x))
        new_z = np.zeros(len(new_x))

        # Z軸については古い値を二乗しておく。
        if self.__err:
            old_z = old_z ** 2

        # バンチング
        for i in range(len(old_x)):
            # 古いX軸の区切りの最小値と最大値を取得(old_mid_x的な)
            # 最小値
            if i == 0:
                x_min = (3*old_x[0]-old_x[1])/2
            else:
                x_min = (old_x[i-1]+old_x[i])/2

            # 最大値
            if i == len(old_x)-1:
                x_max = (3*old_x[-1]-old_x[-2])/2
            else:
                x_max = (old_x[i+1]+old_x[i])/2

            # メイン処理部分
            if x_max > mid_x[0] and x_min < mid_x[-1]:  # 論外な範囲の場合を除外
                min = -1
                max = -1
                # 古い区切りが入る場所を調べる
                for j in range(len(mid_x)-1):
                    if mid_x[j] <= x_min and x_min < mid_x[j+1]:
                        min = j
                    if mid_x[j] <= x_max and x_max < mid_x[j+1]:
                        max = j

                if min == -1 and max == -1:  # ありえない
                    print('Bug')
                    sys.exit()
                elif min == -1:  # 最小値のみ新しい領域をオーバーしている場合
                    for m in range(0, max, 1):
                        new_y[m] += old_y[i]*(mid_x[m+1]-mid_x[m])/(x_max-x_min)
                        if self.__err:
                            new_z[m] += old_z[i]*(mid_x[m+1]-mid_x[m])/(x_max-x_min)
                    new_y[max] += old_y[i]*(x_max-mid_x[max])/(x_max-x_min)
                    if self.__err:
                        new_z[max] += old_z[i]*(x_max-mid_x[max])/(x_max-x_min)
                elif max == -1:  # 最大値のみ新しい領域をオーバーしている場合
                    new_y[min] += old_y[i]*(mid_x[min+1]-x_min)/(x_max-x_min)
                    if self.__err:
                        new_z[min] += old_z[i]*(mid_x[min+1]-x_min)/(x_max-x_min)
                    for m in range(min+1, len(mid_x)-1, 1):
                        new_y[m] += old_y[i]*(mid_x[m+1]-mid_x[m])/(x_max-x_min)
                        if self.__err:
                            new_z[m] += old_z[i]*(mid_x[m+1]-mid_x[m])/(x_max-x_min)
                elif min == max:  # 新しい領域の１つの範囲に収まってしまっている場合
                    new_y[min] += old_y[i]
                    if self.__err:
                        new_z[min] += old_z[i]
                else:  # 新しい領域の複数の範囲にまたがる場合
                    new_y[min] += old_y[i]*(mid_x[min+1]-x_min)/(x_max-x_min)
                    if self.__err:
                        new_z[min] += old_z[i]*(mid_x[min+1]-x_min)/(x_max-x_min)
                    for m in range(min+1, max, 1):
                        new_y[m] += old_y[i]*(mid_x[m+1]-mid_x[m])/(x_max-x_min)
                        if self.__err:
                            new_z[m] += old_z[i]*(mid_x[m+1]-mid_x[m])/(x_max-x_min)
                    new_y[max] += old_y[i]*(x_max-mid_x[max])/(x_max-x_min)
                    if self.__err:
                        new_z[max] += old_z[i]*(x_max-mid_x[max])/(x_max-x_min)

        # insert data
        if self.__err:
            self.data = np.array([new_x, new_y, np.sqrt(new_z)]).T
        else:
            self.data = np.array([new_x, new_y]).T


class MPAList:
    """Analyze 2D List from MPA-3"""
    __filename = ''

    def __init__(self, filename):
        """filename : list file generated with anampa."""
        self.__filename = filename

    def extractTOF(self, ADC=True):
        """Extract Time of Flight (not gated)"""
        if os.path.exists(self.__filename):
            pass
        else:
            print('File %s does not exist' % self.__filename)
            sys.exit()
        outname = os.path.basename(self.__filename).rsplit('.')[0]+'.tof'
        o = open(outname, 'w')
        f = open(self.__filename, 'rb')
        sp = [0]*2048
        while True:
            bytes = f.read(4)
            if bytes == b'':  # バイナリファイルの終端に達したことを示す
                break
            else:
                if ADC:
                    tof, ph = struct.unpack('hh', bytes)
                else:
                    ph, tof = struct.unpack('hh', bytes)
                sp[tof] += 1
        for i in sp:
            o.write('%d\n' % i)
        o.close
        f.close

    def extractPH(self, ADC=True):
        """Extract PulseHeight (not gated)"""
        if os.path.exists(self.__filename):
            pass
        else:
            print('File %s does not exist' % self.__filename)
            sys.exit()
        outname = os.path.basename(self.__filename).rsplit('.')[0]+'.ph'
        o = open(outname, 'w')
        f = open(self.__filename, 'rb')
        sp = [0]*2048
        while True:
            bytes = f.read(4)
            if bytes == b'':  # バイナリファイルの終端に達したことを示す
                break
            else:
                if ADC:
                    tof, ph = struct.unpack('hh', bytes)
                else:
                    ph, tof = struct.unpack('hh', bytes)
                sp[ph] += 1
        for i in sp:
            o.write('%d\n' % i)
        o.close
        f.close

    def gatePH(self, tof_low, tof_high, filename=None, ADC=True):
        """Extract PulseHeight
        1st arg: the lower gate
        2nd arg: the higher gate
        """
        if os.path.exists(self.__filename):
            pass
        else:
            print('File %s does not exist' % self.__filename)
            sys.exit()
        if filename is None:
            outname = (os.path.basename(self.__filename).rsplit('.')[0] + '_' + str(tof_low) + '-' + str(tof_high) + '.ph')
        else:
            outname = filename
        o = open(outname, 'w')
        f = open(self.__filename, 'rb')
        sp = [0]*2048
        while True:
            bytes = f.read(4)
            if bytes == b'':  # バイナリファイルの終端に達したことを示す
                break
            else:
                if ADC:
                    tof, ph = struct.unpack('hh', bytes)
                else:
                    ph, tof = struct.unpack('hh', bytes)
                if tof >= tof_low and tof <= tof_high:  # 境界の値に注意
                    sp[ph] += 1
        for i in sp:
            o.write('%d\n' % i)
        o.close
        f.close

    def ext2d(self):
        """Draw 2d graph"""
        if os.path.exists(self.__filename):
            pass
        else:
            print('File %s does not exist' % self.__filename)
            sys.exit()
        outname = (os.path.basename(self.__filename).rsplit('.')[0] + '_2d.png')
        print("Produced Filename: %s" % outname)
        ct = np.zeros((2048, 2048))
        f = open(self.__filename, 'rb')
        while True:
            bytes = f.read(4)
            if bytes == b'':  # バイナリファイルの終端に達したことを示す
                break
            else:
                tof, ph = struct.unpack('hh', bytes)
                ct[tof-1][ph-1] += 1
        f.close
        # そのままプロットすると、ピーク位置などの色が強く出過ぎてしまう
        # 見た目を良くするためにカウントの対数を取る。
        # カウントが無いチャンネルは-1にする
        # そこそこきれいになる
        for i in range(len(ct)):
            for j in range(len(ct[0])):
                if ct[i][j] == 0:
                    ct[i][j] = -1
                else:
                    ct[i][j] = np.log10(ct[i][j])
        plt.title("%s" % self.__filename)
        plt.imshow(ct)
        plt.savefig(outname)


def ch2keV(ch_width, time_cal, flight_path):
    """calculates neutron Energy from channel width for 15-100keV experiment
    1st arg : channel width between (p, g) peak to (p,n) channel
    2nd arg : time_cal [ns/ch]  (may be about 0.23)
    3rd arg : flight_path[m]    (may be about 0.3)
    """
    c = 0.299792458          # [m/ns]
    time_of_flight = time_cal*ch_width + flight_path/c
    neutron_energy = (72.3 * flight_path / time_of_flight) ** 2
    return neutron_energy * 1000


def ch2keV550(ch_width, time_cal, flight_path):
    """calculates neutron Energy from channel width for 550keV experiment
    1st arg : channel width between (p, g) peak to (p,n) channel
    2nd arg : time_cal [ns/ch]  (may be about 0.23)
    3rd arg : flight_path[m]    (may be about 4.6)
    """
    c = 0.299792458          # [m/ns]
    time_of_flight = time_cal*ch_width + flight_path/c + 500
    neutron_energy = (72.3 * flight_path / time_of_flight) ** 2
    return neutron_energy * 1000


def determineAverageAngle(sampleDiameter, flight_path):
    """ Returns the angle which the sample takes average neutron
    1st arg : sampleDiameter
    2nd arg : flight_path
    Both args must be same unit."""
    return 180*np.arctan(sampleDiameter/3/flight_path)/np.pi


def extractADC1(ifile, ofile):
    """ Extract ADC1 data from MPA file.
    1st arg : input file
    2nd arg : output file
    """
    try:
        outfile = open(ofile, 'w')
        adc1flag = False
        for l in open(ifile):
            if l == '[DATA1,2048 ]\n' or 'CDAT' in l:
                adc1flag = False
            if adc1flag:
                outfile.write(l)
            if l == '[DATA0,2048 ]\n':
                adc1flag = True
                outfile.close
    except:
        print("Cannot open %s" % ifile)


def extractADC2(ifile, ofile):
    """ Extract ADC2 data from MPA file.
    1st arg : input file
    2nd arg : output file
    """
    try:
        outfile = open(ofile, 'w')
        adc2flag = False
        for l in open(ifile):
            if 'CDAT' in l:
                adc2flag = False
            if adc2flag:
                outfile.write(l)
            if l == '[DATA1,2048 ]\n':
                adc2flag = True
                outfile.close
    except:
        print("Cannot open %s" % ifile)


def gauss_func(p, x, y):
    """
    1st arg: parameters
    x, y: expelimental data
    return : residual
    """
    # fx = p[0] * np.exp(-(x-p[1])**2/2.0/p[2]**2)+p[3]*x+p[4]
    fx = p[0] * np.exp(-(x-p[1])**2/2.0/p[2]**2)+p[3]
    residual = y - fx
    return residual


def plotme(gnuplottext, filename=None):
    '''give string to gnuplot'''
    if filename:
        fout = open(filename, "w")
        fout.write(gnuplottext)
        fout.close()
        gnusubprocess = subprocess.Popen(['gnuplot', filename])
        gnusubprocess.wait()
    else:
        gnusubprocess = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
        gnusubprocess.communicate(bytes(gnuplottext, 'UTF-8'))
        gnusubprocess.wait()


def fit_gauss(filename, ini_param, plot=True):
    # if len(ini_param) != 5:
    if len(ini_param) != 4:
        print("Invalid number of initial parameters.")
        return
    print("-----")
    print(filename)
    signal = []
    for line in open(filename):
        signal.append(float(line))

    ini_ch = ini_param[1]
    ini_sig = ini_param[2]
    Px = np.array(range(int(ini_ch-5*ini_sig), int(ini_ch+5*ini_sig)))
    Py = np.array(signal[int(ini_ch-5*ini_sig):int(ini_ch+5*ini_sig)])

    result, error = scipy.optimize.leastsq(gauss_func, ini_param, args=(Px, Py))
    # amp, ch, sig, bg1, bg2 = result
    amp, ch, sig, bg1 = result
    print("amp: %12.4f" % amp)
    print("ch : %12.4f" % ch)
    print("sig: %12.4f" % sig)
    print("bg1: %12.4f" % bg1)
    # print("bg2: %12.4f" % bg2)

    gnutext = """
#!/usr/local/bin/gnuplot
set terminal aqua
set xrange[0:1000]
set sample 1000
fitfunc(x) = (%12.4f)*exp(-(x-(%12.4f))**2/2.0/(%12.4f)**2)+(%12.4f)
""" % (amp, ch, sig, bg1) + """
set title "%s"
""" % (filename) + """
plot "%s" with histeps
replot fitfunc(x)
""" % filename
    if plot:
        plotme(gnutext, filename=filename + ".plt")
    return ch


def rebin(filename, min, max, step, outname=None):
    '''
    Useful function to rebin.
    1st arg: input filename (must be 2 or 3 column data)
    2nd arg: minimum value of new x-region
    3rd arg: maximum value of new x-region
    4th arg: step of new x-region
    outname: Optional. Set your own output-filename.

    Example :
    rebin("background.dat", 0, 10000, 50)
    rebin("neutron_spectrum.dat", 0, 100, 1, outname="rebined_neutron_spectrum.dat")
    '''
    b = Rebin(filename)
    b.rebin(min, max, step)
    if outname:
        b.save(outname)
    else:
        root, ext = os.splitext(filename)
        outname = root + "_rebinned" + ext
        b.save(outname)
