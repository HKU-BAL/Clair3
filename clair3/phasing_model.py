import torch
from torch import nn
from torch.nn import functional as F

import shared.param_f as param


class PhasingModel(nn.Module):
    """Predict per-read haplotype assignment from a phasing matrix.

    The phasing matrix M[R, V] encodes, for each read and each nearby
    heterozygous SNP, whether the read supports the reference allele (1),
    the alternate allele (2), is ambiguous (3), or has no coverage (0).

    The model learns to partition reads into two haplotype groups.

    Architecture:
        1. Embedding: allele categories (4) -> d_emb
        2. Per-read MLP over variant dimension: [R, V*d_emb] -> [R, d_hidden]
        3. Cross-read self-attention (TransformerEncoder): [R, d_hidden] -> [R, d_hidden]
        4. Output head: Linear -> sigmoid -> H in [0,1]^R
    """

    def __init__(self, max_variants=None, d_emb=16, d_hidden=64,
                 n_attn_heads=4, n_attn_layers=2, matrix_depth=None, dropout=0.1):
        super().__init__()
        if max_variants is None:
            max_variants = param.MAX_NEARBY_HETE_SNPS
        if matrix_depth is None:
            matrix_depth = param.matrix_depth_dict.get('ont', 89)

        self.max_variants = max_variants
        self.matrix_depth = matrix_depth
        self.d_hidden = d_hidden

        # 4 allele categories: 0=no_coverage (padding), 1=ref, 2=alt, 3=ambiguous
        self.allele_embed = nn.Embedding(4, d_emb, padding_idx=0)

        self.variant_mlp = nn.Sequential(
            nn.Linear(max_variants * d_emb, d_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_hidden, d_hidden),
        )

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_hidden,
            nhead=n_attn_heads,
            dim_feedforward=d_hidden * 2,
            dropout=dropout,
            batch_first=True,
        )
        self.cross_read_attn = nn.TransformerEncoder(
            encoder_layer, num_layers=n_attn_layers
        )

        self.output_head = nn.Linear(d_hidden, 1)

    def forward(self, M, read_mask=None):
        """
        Args:
            M: [B, R, V] long tensor, allele states (0-3).
            read_mask: [B, R] bool tensor, True for padded/absent reads
                       that should be ignored in attention.

        Returns:
            H: [B, R] float tensor in [0, 1].
               H ~ 0 means HP1, H ~ 1 means HP2, H ~ 0.5 means uncertain.
        """
        B, R, V = M.shape
        x = self.allele_embed(M)          # [B, R, V, d_emb]
        x = x.view(B, R, -1)             # [B, R, V*d_emb]

        # Pad or truncate variant dimension to match expected size
        expected = self.max_variants * self.allele_embed.embedding_dim
        if x.shape[2] < expected:
            x = F.pad(x, (0, expected - x.shape[2]))
        elif x.shape[2] > expected:
            x = x[:, :, :expected]

        x = self.variant_mlp(x)           # [B, R, d_hidden]

        if read_mask is not None:
            x = self.cross_read_attn(x, src_key_padding_mask=read_mask)
        else:
            x = self.cross_read_attn(x)

        h = torch.sigmoid(self.output_head(x).squeeze(-1))  # [B, R]
        return h


class PhasingLoss(nn.Module):
    """Symmetry-aware binary cross-entropy for haplotype prediction.

    Since the assignment of which haplotype is "1" vs "2" is arbitrary,
    we compute BCE for both label orientations and take the minimum.
    Only phased reads (hp_labels > 0) contribute to the loss.
    """

    def forward(self, H, hp_labels, read_mask=None):
        """
        Args:
            H: [B, R] predicted haplotype probability (0=HP1, 1=HP2).
            hp_labels: [B, R] ground truth (0=unphased, 1=HP1, 2=HP2).
            read_mask: [B, R] bool, True for padding reads to ignore.

        Returns:
            Scalar loss.
        """
        # Only compute loss on phased reads
        phased = hp_labels > 0
        if read_mask is not None:
            phased = phased & ~read_mask

        if phased.sum() == 0:
            return torch.tensor(0.0, device=H.device, requires_grad=True)

        H_clamped = torch.clamp(H, min=1e-7, max=1 - 1e-7)

        # Forward orientation: HP1->0, HP2->1
        target_fwd = (hp_labels == 2).float()
        # Reversed orientation: HP1->1, HP2->0
        target_rev = (hp_labels == 1).float()

        bce_fwd = F.binary_cross_entropy(H_clamped, target_fwd, reduction='none')
        bce_rev = F.binary_cross_entropy(H_clamped, target_rev, reduction='none')

        loss_fwd = (bce_fwd * phased.float()).sum() / phased.float().sum()
        loss_rev = (bce_rev * phased.float()).sum() / phased.float().sum()

        return torch.min(loss_fwd, loss_rev)


def inject_hp_channel(tensor, H):
    """Differentiable HP channel injection.

    Maps predicted haplotype probability H in [0,1] to HP channel values
    matching Clair3's encoding: HP1=30, unphased=60, HP2=90.

    The mapping is: hp_value = 30 + H * 60, which is differentiable w.r.t. H.

    Args:
        tensor: [B, R, P, C] float tensor (FA tensor, already float).
                Channel index 7 is the phasing_info channel.
        H: [B, R] float tensor in [0, 1].

    Returns:
        tensor with channel 7 replaced by predicted HP values.
        Only positions where the read exists (mapping_quality channel != 0)
        are modified.
    """
    hp_values = 30.0 + H * 60.0             # [B, R], range [30, 90]
    hp_expanded = hp_values.unsqueeze(2)     # [B, R, 1]
    hp_expanded = hp_expanded.expand(-1, -1, tensor.shape[2])  # [B, R, P]

    # Mask: only set HP where read actually exists (mapping_quality != 0)
    read_exists = (tensor[:, :, :, 2] != 0).float()  # [B, R, P]

    tensor = tensor.clone()
    tensor[:, :, :, 7] = hp_expanded * read_exists
    return tensor


def clear_hp_channel(tensor):
    """Set HP channel to unphased (60) for all positions where read exists.

    Used in frozen/joint modes to remove WhatsHap HP information before
    injecting neural phasing predictions.

    Args:
        tensor: [B, R, P, C] float tensor.

    Returns:
        tensor with channel 7 set to 60 (unphased) where read exists, 0 elsewhere.
    """
    tensor = tensor.clone()
    read_exists = (tensor[:, :, :, 2] != 0).float()
    tensor[:, :, :, 7] = 60.0 * read_exists
    return tensor
