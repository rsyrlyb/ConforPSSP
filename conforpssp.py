from tensorflow.keras.models import load_model
import scripts as mf
import numpy as np

MAX_TEXT_LEN = 202
ltok =4


class PSSP:
    def __init__(self, n):
        self.n = n
        #self.model_name = 'a_' + n + '_20_202_model'
        self.model_name = "model_" + n

        self.model = load_model("models/" + self.model_name)
        frd = mf.readjson('db/freq_table.json')
        aad = mf.readjson('db/aa_tokens_dict.json')
        aas = mf.readjson('db/ss_tokens_dict.json')
        self.frd = frd

        self.word_dict_inp = mf.modify_dict(aad)
        self.word_dict_out = mf.modify_dict(aas)
        self.word_dict_outr = {self.word_dict_out[x]: x for x in self.word_dict_out}
        self.word_dict_inpr = {self.word_dict_inp[x]: x for x in self.word_dict_inp}

    def generate(self, aa_seq, aa_file, output_dir, num_iter=1, k_prob=10):

        output = []
        lt = (MAX_TEXT_LEN - 2) * ltok
        self.k_prob = k_prob
        aa_seqs = {'str': aa_seq}
        if aa_file is not None:
            aa_list = [a.split('\n') for a in mf.readfile(aa_file,'t').split('>')]
            aa_seqs = {a[0]: a[1] for a in aa_list if len(a[:2]) == 2}

        str_list = [a for a in aa_seqs]
        str_list_fr = mf.make_ranges(str_list, min(len(str_list), 20))
        for sx, str_list in enumerate(str_list_fr):
            print(sx)
            inp = []
            #unq_res = {seq + '_' + str(ax) + '_': [] for ax, seq in enumerate(str_list)}
            unq_res = {seq + '_': [] for ax, seq in enumerate(str_list)}
      
            for ax, seq in enumerate(str_list):
                aa_seq = aa_seqs[seq][:lt]
                x_pos = []
                if 'X' in aa_seq:
                    x_pos = [ax for ax in range(len(aa_seq)) if aa_seq[ax] == 'X']
                for it in range(num_iter):
                    i_name = '_'.join([seq, str(it)])
                    #i_name = '_'.join([seq, str(ax), str(it)])
                    inp.append([i_name, x_pos, self.input_processing(aa_seq)])

       
            i_names = [a[0] for a in inp]
            x_pos = [a[1] for a in inp]
            inp = [a[2] for a in inp]
            inp = [np.vstack([inp[y][x] for y in range(len(inp))]) for x in range(len(inp[0]))]
            model_output = self.model.predict(inp[:], verbose=None)
            inp_aa = mf.make_ranges(inp[0], 4)
            inp_fr = mf.make_ranges(inp[1], 4)
            model_output = mf.make_ranges(model_output, 4)
            for x in range(len(model_output)):
                inp_pss, pss = self.output_processing(model_output[x], inp_aa[x], x_pos[x])
                if pss not in unq_res['_'.join(i_names[x].split('_')[:-1]) + '_']:
                    fr_pss = self.fr_processing(inp_fr[x][ltok-1])
                    #if pss not in [v[1] for v in output]:
                    output.append([i_names[x], pss, fr_pss])
                    #unq_res['_'.join(i_names[x].split('_')[:-1]) + '_'].append(pss)
            #output.append(output_ax)
        output = '\n'.join(['>' + a[0] + '\n' + a[1] for a in output]) #+ '\n' + a[2] 
        #output = '\n'.join(['>' + a[0] + '\n' + a[1] + '\n' + a[2] for a in output])


        if aa_file is None:
            output_file_name = output_dir + 'pssp_' + self.n + '.fasta'
        else:
            aa_file = aa_file.split('/')[-1]
            output_file_name = output_dir + aa_file.replace('.fasta', '_pssp_' + self.model_name + '.fasta')
        print(output_file_name)
        mf.writefile(output_file_name, output)
        return output

    def fr_processing(self, fr_seq):
  
        end_idx = [ax for ax in range(len(fr_seq)) if fr_seq[ax] == 2][0]
        fr_seq = fr_seq[1:end_idx]
        gap_ind = [ax for ax in range(len(fr_seq)) if fr_seq[ax] == mf.word_dict['[GAP]']]
        fr_seq = [a if ax not in gap_ind else mf.word_dict['[GAP]'] for ax, a in enumerate(fr_seq)]
        fr_seq = ''.join([self.word_dict_outr[a] for a in fr_seq]).replace('[GAP]', '!')
        return fr_seq

    def input_processing(self, aa_seq):
        te = 'X' * (ltok - 1)
        #aa_seq = aa_seq.replace('J', '!')
        aa_seq_red = 'J'.join([te + i + te for i in mf.rC_red_aa(aa_seq).split('!')])
        input_batch = mf.generate_inputs(aa_seq_red, self.frd, self.word_dict_inp, self.word_dict_out,
                                         self.word_dict_outr, MAX_TEXT_LEN, self.k_prob)
        return input_batch

    def output_processing(self, tar_pred, tokens_a, x_pos):

        end_idx = [[ax for ax in range(len(tokens_a[x])) if tokens_a[x][ax] == 2][0] for x in range(ltok)]
        tar_pred = [tar_pred[x][:end_idx[x] - 1] for x in range(ltok)]
        max_v = np.array([tar_pred[k].max(axis=1).mean() for k in range(ltok)]).argmax()
        tar_pred = tar_pred[max_v].argmax(axis=1)

        tokens_a = [tokens_a[x][1:end_idx[x]] for x in range(ltok)]
        tokens_a = tokens_a[max_v]
        gap_ind = [ax for ax in range(len(tokens_a)) if tokens_a[ax] == mf.word_dict['[GAP]']]

        tar_pred = [a if ax not in gap_ind else mf.word_dict['[GAP]'] for ax, a in enumerate(tar_pred)]
        tar_pred = ''.join([self.word_dict_outr[a] for a in tar_pred]).replace('[GAP]', '!')

        tokens_a =  [a if ax not in gap_ind else mf.word_dict['[GAP]'] for ax, a in enumerate(tokens_a)]
        inp = ''.join([self.word_dict_inpr[a] for a in tokens_a]).replace('[GAP]', '!')

        kx = [ax for ax in range(len(inp)) if (inp[ax] == 'X')]
        tar_pred = ''.join([a for ax, a in enumerate(tar_pred) if ax not in kx])
        inp = ''.join([a for ax, a in enumerate(inp) if ax not in kx])
        if x_pos != []:
            for pos in x_pos:
                tar_pred = tar_pred[:pos] + 'X' + tar_pred[pos:]
                inp = inp[:pos] + 'X' + inp[pos:]
        return inp, tar_pred




