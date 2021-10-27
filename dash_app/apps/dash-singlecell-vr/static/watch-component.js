/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('watch', {
    schema: {
        targetId: {type: "string"}
    },
    init: function() {
        this.tock = AFRAME.utils.throttleTick(this.tock, 500, this);
        this.vector = new THREE.Vector3();
        this.rig = document.getElementById('rig').object3D;
        this.campos = rig.position;
    },
    tock: function() {
        if (this.campos !== this.rig.position) {
            const target = this.el.sceneEl.camera;

            target.updateMatrixWorld();
            this.vector.setFromMatrixPosition(target.matrixWorld);

            if (this.el.object3D.parent) {
                this.el.object3D.parent.updateMatrixWorld();
                this.el.object3D.parent.worldToLocal(this.vector);
            }

            this.el.object3D.lookAt(this.vector);
        }
    }
});