'''
    References
    ----------
    @misc{https://doi.org/10.48550/arxiv.1706.03762,
        doi = {10.48550/ARXIV.1706.03762},
        url = {https://arxiv.org/abs/1706.03762},
        author = {Vaswani, Ashish and Shazeer, Noam and Parmar, Niki and Uszkoreit, Jakob and Jones, Llion and Gomez, Aidan N. and Kaiser, Lukasz and Polosukhin, Illia},
        keywords = {Computation and Language (cs.CL), Machine Learning (cs.LG), FOS: Computer and information sciences, FOS: Computer and information sciences},
        title = {Attention Is All You Need},
        publisher = {arXiv},
        year = {2017},
        copyright = {arXiv.org perpetual, non-exclusive license}
    }

    @inproceedings{wolf-etal-2020-transformers,
        title = "Transformers: State-of-the-Art Natural Language Processing",
        author = "Thomas Wolf and Lysandre Debut and Victor Sanh and Julien Chaumond and Clement Delangue and Anthony Moi and Pierric Cistac and Tim Rault and Rémi Louf and Morgan Funtowicz and Joe Davison and Sam Shleifer and Patrick von Platen and Clara Ma and Yacine Jernite and Julien Plu and Canwen Xu and Teven Le Scao and Sylvain Gugger and Mariama Drame and Quentin Lhoest and Alexander M. Rush",
        booktitle = "Proceedings of the 2020 Conference on Empirical Methods in Natural Language Processing: System Demonstrations",
        month = oct,
        year = "2020",
        address = "Online",
        publisher = "Association for Computational Linguistics",
        url = "https://www.aclweb.org/anthology/2020.emnlp-demos.6",
        pages = "38--45"
    }

'''

# Code:   jiaqili@zju.edu
#         https://github.com/beiweixiaoxu/transformerencoder/blob/master/TransformerEncoder.py
# Note:   modified for onehot sequence input

import logging
import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from .BasicModule import BasicModule

class ScaledDotProductAttention(nn.Module):
    def __init__(self, d_k):
        super(ScaledDotProductAttention, self).__init__()
        self.d_k = d_k

    def forward(self, q, k, v, attn_mask):
        # |q| : (batch_size, n_heads, q_len, d_k), |k| : (batch_size, n_heads, k_len, d_k), |v| : (batch_size, n_heads, v_len, d_v)
        # |attn_mask| : (batch_size, n_heads, seq_len(=q_len), seq_len(=k_len))
        attn_score = torch.matmul(q, k.transpose(-1, -2)) / np.sqrt(self.d_k)
        attn_score.masked_fill_(attn_mask, -1e9)
        # |attn_score| : (batch_size, n_heads, q_len, k_len)
        attn_weights = nn.Softmax(dim=-1)(attn_score)
        # |attn_weights| : (batch_size, n_heads, q_len, k_len)
        output = torch.matmul(attn_weights, v)
        # |output| : (batch_size, n_heads, q_len, d_v)
        return output, attn_weights


class MultiHeadAttention(nn.Module):
    def __init__(self, d_model, n_heads):
        super(MultiHeadAttention, self).__init__()
        self.n_heads = n_heads
        self.d_k = self.d_v = d_model // n_heads // 4
        self.WQ = nn.Linear(d_model, d_model // 4)
        self.WK = nn.Linear(d_model, d_model // 4)
        self.WV = nn.Linear(d_model, d_model // 4)
        self.scaled_dot_product_attn = ScaledDotProductAttention(self.d_k)
        self.linear = nn.Linear(n_heads * self.d_v, d_model)

    def forward(self, Q, K, V, attn_mask):
        # |Q| : (batch_size, q_len, d_model), |K| : (batch_size, k_len, d_model), |V| : (batch_size, v_len, d_model)
        # |attn_mask| : (batch_size, seq_len(=q_len), seq_len(=k_len))
        batch_size = Q.size(0)
        q_heads = self.WQ(Q).view(batch_size, -1, self.n_heads, self.d_k).transpose(1, 2)
        k_heads = self.WK(K).view(batch_size, -1, self.n_heads, self.d_k).transpose(1, 2)
        v_heads = self.WV(V).view(batch_size, -1, self.n_heads, self.d_v).transpose(1, 2)
        # |q_heads| : (batch_size, n_heads, q_len, d_k), |k_heads| : (batch_size, n_heads, k_len, d_k), |v_heads| : (batch_size, n_heads, v_len, d_v)
        attn_mask = attn_mask.unsqueeze(1).repeat(1, self.n_heads, 1, 1)
        # |attn_mask| : (batch_size, n_heads, seq_len(=q_len), seq_len(=k_len))
        attn, attn_weights = self.scaled_dot_product_attn(q_heads, k_heads, v_heads, attn_mask)
        # |attn| : (batch_size, n_heads, q_len, d_v)
        # |attn_weights| : (batch_size, n_heads, q_len, k_len)
        attn = attn.transpose(1, 2).contiguous().view(batch_size, -1, self.n_heads * self.d_v)
        # |attn| : (batch_size, q_len, n_heads * d_v)
        output = self.linear(attn)
        # |output| : (batch_size, q_len, d_model)
        return output, attn_weights


class PositionWiseFeedForwardNetwork(nn.Module):
    def __init__(self, d_model, d_ff):
        super(PositionWiseFeedForwardNetwork, self).__init__()
        self.linear1 = nn.Linear(d_model, d_ff)
        self.linear2 = nn.Linear(d_ff, d_model)
        self.relu = nn.ReLU()

    def forward(self, inputs):
        # |inputs| : (batch_size, seq_len, d_model)
        output = self.relu(self.linear1(inputs))
        # |output| : (batch_size, seq_len, d_ff)
        output = self.linear2(output)
        # |output| : (batch_size, seq_len, d_model)
        return output


class EncoderLayer(nn.Module):
    def __init__(self, d_model, n_heads, p_drop, d_ff):
        super(EncoderLayer, self).__init__()
        self.mha = MultiHeadAttention(d_model, n_heads)
        self.dropout1 = nn.Dropout(p_drop)
        self.layernorm1 = nn.LayerNorm(d_model, eps=1e-6)
        self.ffn = PositionWiseFeedForwardNetwork(d_model, d_ff)
        self.dropout2 = nn.Dropout(p_drop)
        self.layernorm2 = nn.LayerNorm(d_model, eps=1e-6)

    def forward(self, inputs, attn_mask):
        # |inputs| : (batch_size, seq_len, d_model)
        # |attn_mask| : (batch_size, seq_len, seq_len)
        attn_outputs, attn_weights = self.mha(inputs, inputs, inputs, attn_mask)
        attn_outputs = self.dropout1(attn_outputs)
        attn_outputs = self.layernorm1(inputs + attn_outputs)
        # |attn_outputs| : (batch_size, seq_len(=q_len), d_model)
        # |attn_weights| : (batch_size, n_heads, q_len, k_len)
        ffn_outputs = self.ffn(attn_outputs)
        ffn_outputs = self.dropout2(ffn_outputs)
        ffn_outputs = self.layernorm2(attn_outputs + ffn_outputs)
        # |ffn_outputs| : (batch_size, seq_len, d_model)
        return ffn_outputs, attn_weights


class TransformerEncoder(BasicModule):
    """Transformer Encoder in NvTK. 
    TransformerEncoder is a stack of MultHeadAttention encoder layers.

    Args:
        vocab_size (int)    : vocabulary size (vocabulary: collection mapping token to numerical identifiers)
        seq_len    (int)    : input sequence length
        d_model    (int)    : number of expected features in the input
        n_layers   (int)    : number of sub-encoder-layers in the encoder
        n_heads    (int)    : number of heads in the multiheadattention models
        p_drop     (float)  : dropout value
        d_ff       (int)    : dimension of the feedforward network model
        pad_id     (int)    : pad token id

    Examples:
    >>> encoder = TransformerEncoder(vocab_size=1000, seq_len=512)
    >>> inp = torch.arange(512).repeat(2, )
    >>> encoder(inp)
    
    """
    def __init__(self, seq_len, vocab_size=4, d_model=512, n_layers=6, n_heads=8, p_drop=0.1, d_ff=2048, pad_id=torch.zeros(4), embedding=None, embedding_weight=None, fix_embedding=False, max_atten=True):
        super(TransformerEncoder, self).__init__()
        # word embedding
        self.pad_id = pad_id
        # embedding layers
        if embedding is not None: 
            self.embedding = embedding
        else:
            self.embedding = nn.Embedding(vocab_size, d_model) 
        if embedding_weight: 
            self.embedding.weight = nn.Parameter(embedding_weight) # 如果有训练好的emb输入
            # 也可以写成 self.emb = nn.from_pretrained(embedding_weight)
        # self.embedding.weight.requires_grad = False if fix_embedding else True # 设定emb是否随训练更新
        
        # pos_embedding
        self.sinusoid_table = self.get_sinusoid_table(seq_len + 1, d_model)  # (seq_len+机器学习, d_model)  # ?开头?
        self.pos_embedding = nn.Embedding.from_pretrained(self.sinusoid_table, freeze=True)

        # EncoderLayer
        self.layers = nn.ModuleList([EncoderLayer(d_model, n_heads, p_drop, d_ff) for _ in range(n_layers)])
        
        # attention_weight
        self.attention_weights = None

        self.max_atten = max_atten

        # layers to classify
        # self.linear = nn.Linear(d_model, 2)  # 不需要加linear层
        # self.log_softmax = F.log_softmax

    def forward(self, inputs):
        # |inputs| : (batch_size, channel=4, seq_len)
        logging.debug(inputs.shape)

        inputs = self.embedding(inputs) # (batch_size, d_model, seq_len)
        outputs = inputs.transpose(1, -1) # (batch_size, seq_len, d_model)
        logging.debug(outputs.shape)
        logging.debug(outputs)

        positions = torch.arange(inputs.size(-1), device=inputs.device, dtype=torch.long).repeat(inputs.size(0), 1) + 1
        # |positions| : (batch_size * seq_len)
        logging.debug(positions.shape)        
        
        position_pad_mask = inputs.sum(1) == 0 #inputs.eq(self.pad_id)  # position_pad_mask:pad过的位置为1，其它位置为0
        positions.masked_fill_(position_pad_mask, 0)  # positions:pad过的位置为0，其它位置按原来的标号
        # |positions| : (batch_size, seq_len)
        logging.debug(positions.shape)
        logging.debug(positions)
        
        outputs = outputs + self.pos_embedding(positions)
        # |outputs| : (batch_size, seq_len, d_model)
        logging.debug(outputs.shape)
        logging.debug(outputs)
        
        attn_pad_mask = self.get_attention_padding_mask(inputs, inputs, self.pad_id)
        # |attn_pad_mask| : (batch_size, seq_len, seq_len)
        logging.debug(attn_pad_mask.shape)
        logging.debug(attn_pad_mask)
        
        attention_weights = []
        for layer in self.layers:
            # |outputs| : (batch_size, seq_len, d_model)
            outputs, attn_weights = layer(outputs, attn_pad_mask)
            # |attn_weights| : (batch_size, n_heads, seq_len, seq_len)
            attention_weights.append(attn_weights.data)
        # |attn_weights| : (n_layers, batch_size, n_heads, seq_len, seq_len)
        self.attention_weights = attention_weights # np.array(attention_weights)

        # |outputs| : (batch_size, seq_len, d_model)
        if self.max_atten:
            outputs, _ = torch.max(outputs, dim=1)
        logging.debug(outputs.shape)
        logging.debug(attention_weights[0].shape)
        # |outputs| : (batch_size, d_model)
        # outputs = self.log_softmax(self.linear(outputs),dim=-1)
        # |outputs| : (batch_size, 2)
        return outputs #, attention_weights

    def get_attention_padding_mask(self, q, k, pad_id):
        """Mask Attention Padding.

        Args:
            q   (torch.Tensor) : query tensor
            k   (torch.Tensor) : key tensor
            pad_id  (int)   : pad token id

        Return:
            attn_pad_mask (torch.BoolTensor)  :   Attention Padding Masks
        """
        attn_pad_mask = (k.sum(1) == 0 ).unsqueeze(1).repeat(1, q.size(-1), 1)
        # k.eq(pad_id).unsqueeze(1).repeat(1, q.size(1), 1)
        # |attn_pad_mask| : (batch_size, q_len, k_len)
        return attn_pad_mask

    def get_sinusoid_table(self, seq_len, d_model):
        """Sinusoid Position encoding table in transformer.

        Args:
            seq_len   (int) : sequence length
            d_model   (int) : model dimension

        Return:
            sinusoid_table (torch.FloatTensor)  :   Sinusoid Position encoding table

        """
        def get_angle(pos, i, d_model):
            return pos / np.power(10000, (2 * (i // 2)) / d_model)

        sinusoid_table = np.zeros((seq_len, d_model))
        for pos in range(seq_len):
            for i in range(d_model):
                if i % 2 == 0:
                    sinusoid_table[pos, i] = np.sin(get_angle(pos, i, d_model))
                else:
                    sinusoid_table[pos, i] = np.cos(get_angle(pos, i, d_model))
        return torch.FloatTensor(sinusoid_table)

    def get_attention(self):
        """Get the attention weights of Transformer Encoder

        Return:
            attention_weights (torch.FloatTensor)  :    attention weights 

        """
        # |attn_weights| : list of (n_layers, batch_size, n_heads, seq_len, seq_len)
        return self.attention_weights


# def get_transformer_attention(model, data_loader, device=torch.device("cuda")):
#     attention = []
    
#     model.eval()
#     for data, target in data_loader:
#         data, target = data.to(device), target.to(device)
#         pred = model(data)
#         batch_attention = model.Embedding.get_attention()
#         batch_attention = np.array([atten.cpu().data.numpy() for atten in batch_attention]).swapaxes(0,1)
#         attention.append(batch_attention)

#     attention = np.concatenate(attention, 0) # (size, n_layers, n_heads,seq_len, seq_len)
#     return attention
